// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUTPCNeighboursFinder.cxx
/// \author Sergey Gorbunov, Ivan Kisel, David Rohr

#include "GPUTPCHit.h"
#include "GPUTPCNeighboursFinder.h"
#include "GPUTPCTracker.h"
#include "GPUCommonMath.h"
#include "GPUDefMacros.h"
using namespace GPUCA_NAMESPACE::gpu;

template <>
GPUdii() void GPUTPCNeighboursFinder::Thread<0>(int /*nBlocks*/, int nThreads, int iBlock, int iThread, GPUsharedref() MEM_LOCAL(GPUSharedMemory) & GPUrestrict() ss, processorType& GPUrestrict() trackerX)
{
  //* find neighbours

#ifdef GPUCA_GPUCODE
  for (unsigned int i = iThread; i < sizeof(MEM_PLAIN(GPUTPCRow)) / sizeof(int); i += nThreads) {
    reinterpret_cast<GPUsharedref() int*>(&ss.mRow)[i] = reinterpret_cast<GPUglobalref() int*>(&trackerX.SliceDataRows()[iBlock])[i];
    if (iBlock >= 2 && iBlock < GPUCA_ROW_COUNT - 2) {
      reinterpret_cast<GPUsharedref() int*>(&ss.mRowUp)[i] = reinterpret_cast<GPUglobalref() int*>(&trackerX.SliceDataRows()[iBlock + 2])[i];
      reinterpret_cast<GPUsharedref() int*>(&ss.mRowDown)[i] = reinterpret_cast<GPUglobalref() int*>(&trackerX.SliceDataRows()[iBlock - 2])[i];
    }
  }
  GPUbarrier();
#endif

  HIPGPUsharedref() MEM_LOCAL(GPUSharedMemory) & GPUrestrict() s = (HIPGPUsharedref() MEM_LOCAL(GPUSharedMemory)&)ss;
  HIPGPUconstantref() processorType& GPUrestrict() tracker = (HIPGPUconstantref() processorType&)trackerX;

#ifdef GPUCA_GPUCODE
  HIPGPUsharedref() const MEM_LOCAL(GPUTPCRow) & GPUrestrict() row = s.mRow;
  HIPGPUsharedref() const MEM_LOCAL(GPUTPCRow) & GPUrestrict() rowUp = s.mRowUp;
  HIPGPUsharedref() const MEM_LOCAL(GPUTPCRow) & GPUrestrict() rowDn = s.mRowDown;
#else
  HIPGPUglobalref() const MEM_GLOBAL(GPUTPCRow) & GPUrestrict() row = tracker.mData.mRows[iBlock];
  HIPGPUglobalref() const MEM_GLOBAL(GPUTPCRow) & GPUrestrict() rowUp = tracker.mData.mRows[iBlock + 2];
  HIPGPUglobalref() const MEM_GLOBAL(GPUTPCRow) & GPUrestrict() rowDn = tracker.mData.mRows[iBlock - 2];
#endif

  if (iThread == 0) {
    s.mIRow = iBlock;
    s.mIRowUp = iBlock + 2;
    s.mIRowDn = iBlock - 2;
    if (s.mIRow < GPUCA_ROW_COUNT) {
      s.mNHits = row.mNHits;
      if ((s.mIRow >= 2) && (s.mIRow <= GPUCA_ROW_COUNT - 3)) {
        // the axis perpendicular to the rows
        const float xDn = rowDn.mX;
        const float x = row.mX;
        const float xUp = rowUp.mX;

        // distance of the rows (absolute and relative)
        s.mUpDx = xUp - x;
        s.mDnDx = xDn - x;
        s.mUpTx = xUp / x;
        s.mDnTx = xDn / x;
      }
    }
  }
  GPUbarrier();

  // local copies

  if ((s.mIRow <= 1) || (s.mIRow >= GPUCA_ROW_COUNT - 2) || (rowUp.mNHits < 0) || (rowDn.mNHits < 0)) {
    const int lHitNumberOffset = row.mHitNumberOffset;
    for (int ih = iThread; ih < s.mNHits; ih += nThreads) {
      tracker.mData.mLinkUpData[lHitNumberOffset + ih] = CALINK_INVAL;
      tracker.mData.mLinkDownData[lHitNumberOffset + ih] = CALINK_INVAL;
    }
    return;
  }

  const float chi2Cut = 3.f * 3.f * 4 * (s.mUpDx * s.mUpDx + s.mDnDx * s.mDnDx);
  // float chi2Cut = 3.*3.*(s.mUpDx*s.mUpDx + s.mDnDx*s.mDnDx ); //SG

  const float y0 = row.mGrid.mYMin;
  const float z0 = row.mGrid.mZMin;
  const float stepY = row.mHstepY;
  const float stepZ = row.mHstepZ;

  const int lHitNumberOffset = row.mHitNumberOffset;
  const int lHitNumberOffsetUp = rowUp.mHitNumberOffset;
  const int lHitNumberOffsetDn = rowDn.mHitNumberOffset;
  const unsigned int lFirstHitInBinOffsetUp = rowUp.mFirstHitInBinOffset;
  const unsigned int lFirstHitInBinOffsetDn = rowDn.mFirstHitInBinOffset;
  HIPGPUglobalref() const calink* GPUrestrict() lFirstHitInBin = (HIPGPUglobalref() const calink*)tracker.mData.mFirstHitInBin;
  HIPGPUglobalref() const cahit2* GPUrestrict() pHitData = (HIPGPUglobalref() const cahit2*)tracker.mData.mHitData;
  const float y0Up = rowUp.mGrid.mYMin;
  const float z0Up = rowUp.mGrid.mZMin;
  const float stepYUp = rowUp.mHstepY;
  const float stepZUp = rowUp.mHstepZ;
  const float y0Dn = rowDn.mGrid.mYMin;
  const float z0Dn = rowDn.mGrid.mZMin;
  const float stepYDn = rowDn.mHstepY;
  const float stepZDn = rowDn.mHstepZ;

  const float kAngularMultiplier = tracker.mConstantMem->param.rec.SearchWindowDZDR;
  const float kAreaSizeY = tracker.mConstantMem->param.rec.NeighboursSearchArea;
  const float kAreaSizeZUp = kAngularMultiplier != 0.f ? (s.mUpDx * kAngularMultiplier) : kAreaSizeY;
  const float kAreaSizeZDn = kAngularMultiplier != 0.f ? (-s.mDnDx * kAngularMultiplier) : kAreaSizeY;
  const float kAreaSlopeZUp = kAngularMultiplier != 0.f ? 1.f : s.mUpTx;
  const float kAreaSlopeZDn = kAngularMultiplier != 0.f ? 1.f : s.mDnTx;

#if GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP < GPUCA_MAXN
  calink neighUp[GPUCA_MAXN - GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP];
  float yzUp[GPUCA_MAXN - GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP];
  float yzUp2[GPUCA_MAXN - GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP];
#endif // GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP > 0

  for (int ih = iThread; ih < s.mNHits; ih += nThreads) {

    HIPGPUglobalref() const cahit2& hitData = pHitData[lHitNumberOffset + ih];
    const float y = y0 + (hitData.x) * stepY;
    const float z = z0 + (hitData.y) * stepZ;

    int nNeighUp = 0;
    float minZ, maxZ, minY, maxY;
    int binYmin, binYmax, binZmin, binZmax;
    int nY;

    {
      const float yy = y * s.mUpTx;
      const float zz = z * kAreaSlopeZUp;
      const float dy = kAreaSizeY;
      const float dz = kAreaSizeZUp;
      minZ = zz - dz;
      maxZ = zz + dz;
      minY = yy - dy;
      maxY = yy + dy;
      reinterpret_cast<const HIPTPCROW(GPUTPCRow)&>(rowUp).Grid().GetBin(minY, minZ, &binYmin, &binZmin);
      reinterpret_cast<const HIPTPCROW(GPUTPCRow)&>(rowUp).Grid().GetBin(maxY, maxZ, &binYmax, &binZmax);
      nY = reinterpret_cast<const HIPTPCROW(GPUTPCRow)&>(rowUp).Grid().Ny();
    }

    for (int k1 = binZmin; k1 <= binZmax && (nNeighUp < GPUCA_MAXN); k1++) {
      int iMin = lFirstHitInBin[lFirstHitInBinOffsetUp + k1 * nY + binYmin];
      int iMax = lFirstHitInBin[lFirstHitInBinOffsetUp + k1 * nY + binYmax + 1];
      GPUCA_UNROLL(U(4), U(2))
      for (int i = iMin; i < iMax && (nNeighUp < GPUCA_MAXN); i++) {
        HIPGPUglobalref() const cahit2& hitDataUp = pHitData[lHitNumberOffsetUp + i];
        GPUTPCHit h;
        h.mY = y0Up + (hitDataUp.x) * stepYUp;
        h.mZ = z0Up + (hitDataUp.y) * stepZUp;

        if (h.mY < minY || h.mY > maxY || h.mZ < minZ || h.mZ > maxZ) {
          continue;
        }

#if GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP < GPUCA_MAXN
#if GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP == 0
        if (true) {
#else
        if (nNeighUp >= (int)GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP) {
#endif
          neighUp[nNeighUp - GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP] = (calink)i;
          yzUp[nNeighUp - GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP] = s.mDnDx * (h.Y() - y);
          yzUp2[nNeighUp - GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP] = s.mDnDx * (h.Z() - z);
        } else
#endif
        {
#if GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP > 0
          s.mB[nNeighUp][iThread] = (calink)i;
          s.mA1[nNeighUp][iThread] = s.mDnDx * (h.Y() - y);
          s.mA2[nNeighUp][iThread] = s.mDnDx * (h.Z() - z);
#endif
        }
        nNeighUp++;
      }
    }

#if GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP > 0
    for (int iUp = nNeighUp; iUp < GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP; iUp++) {
      s.mA1[iUp][iThread] = -1.e10f;
      s.mA2[iUp][iThread] = -1.e10f;
      s.mB[iUp][iThread] = (calink)-1;
    }
#endif
    /* 
#if GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP < GPUCA_MAXN
    {
      int N = nNeighUp - GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP;
      int N1 = (N / 4) * 4;
      if (N1 != N) {
        if (N1 + 1 >= N) {
          yzUp[N1 + 1] = -1.e10f;
          yzUp2[N1 + 1] = -1.e10f;
          neighUp[N1 + 1] = (calink)-1;
        }
        if (N1 + 2 >= N) {
          yzUp[N1 + 2] = -1.e10f;
          yzUp2[N1 + 2] = -1.e10f;
          neighUp[N1 + 2] = (calink)-1;
        }
        if (N1 + 3 >= N) {
          yzUp[N1 + 3] = -1.e10f;
          yzUp2[N1 + 3] = -1.e10f;
          neighUp[N1 + 3] = (calink)-1;
        }
      }
    }
#endif
*/

    {
      const float yy = y * s.mDnTx;
      const float zz = z * kAreaSlopeZDn;
      const float dy = kAreaSizeY;
      const float dz = kAreaSizeZDn;
      minZ = zz - dz;
      maxZ = zz + dz;
      minY = yy - dy;
      maxY = yy + dy;
    }
    reinterpret_cast<const HIPTPCROW(GPUTPCRow)&>(rowDn).Grid().GetBin(minY, minZ, &binYmin, &binZmin);
    reinterpret_cast<const HIPTPCROW(GPUTPCRow)&>(rowDn).Grid().GetBin(maxY, maxZ, &binYmax, &binZmax);
    nY = reinterpret_cast<const HIPTPCROW(GPUTPCRow)&>(rowDn).Grid().Ny();

    int linkUp = -1;
    int linkDn = -1;
    int bestUp = -1;
    int bestDn = -1;
    float bestD = 1.e10f;

    for (int k1 = binZmin; k1 <= binZmax; k1++) {
      int iMin = lFirstHitInBin[lFirstHitInBinOffsetDn + k1 * nY + binYmin];
      int iMax = lFirstHitInBin[lFirstHitInBinOffsetDn + k1 * nY + binYmax + 1];
      for (int i = iMin; i < iMax; i++) {
        HIPGPUglobalref() const cahit2& hitDataDn = pHitData[lHitNumberOffsetDn + i];
        GPUTPCHit h;
        h.mY = y0Dn + (hitDataDn.x) * stepYDn;
        h.mZ = z0Dn + (hitDataDn.y) * stepZDn;

        if (h.mY < minY || h.mY > maxY || h.mZ < minZ || h.mZ > maxZ) {
          continue;
        }

        float2 yzdn = CAMath::MakeFloat2(s.mUpDx * (h.Y() - y), s.mUpDx * (h.Z() - z));
/* // original
        GPUCA_UNROLL(U(4), )
        for (int iUp = 0; iUp < nNeighUp; iUp++) {
#if GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP > 0 && GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP < GPUCA_MAXN
          float2 yzup = iUp >= GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP ? CAMath::MakeFloat2(yzUp[iUp - GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP], yzUp2[iUp - GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP]) : CAMath::MakeFloat2(s.mA1[iUp][iThread], s.mA2[iUp][iThread]);
#elif GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP == GPUCA_MAXN
          const float2 yzup = CAMath::MakeFloat2(s.mA1[iUp][iThread], s.mA2[iUp][iThread]);
#else
          const float2 yzup = CAMath::MakeFloat2(yzUp[iUp], yzUp2[iUp]);
#endif
          const float dy = yzdn.x - yzup.x;
          const float dz = yzdn.y - yzup.y;
          const float d = dy * dy + dz * dz;
          if (d < bestD) {
            bestD = d;
            bestDn = i;
            bestUp = iUp;
          }
        }
        */
#if GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP > 0
        constexpr int tmpMaxUp = GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP;
#pragma unroll(tmpMaxUp)
        for (int iUp = 0; iUp < GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP; iUp++) {
          const float dy = yzdn.x - s.mA1[iUp][iThread];
          const float dz = yzdn.y - s.mA2[iUp][iThread];
          const float d = dy * dy + dz * dz;
          if (d < bestD) {
            bestD = d;
            bestDn = i;
            bestUp = iUp;
          }
        }
#endif

#if 0 && GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP < GPUCA_MAXN
        {
          int N = nNeighUp - GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP;
#pragma unroll(4)
          for (int iUp = 0; iUp < N; iUp++) {
            const float dy = yzdn.x - yzUp[iUp];
            const float dz = yzdn.y - yzUp2[iUp];
            const float d = dy * dy + dz * dz;
            if (d < bestD) {
              bestD = d;
              bestDn = i;
              bestUp = GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP + iUp;
            }
          }
        }
#endif

#if 1 && GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP < GPUCA_MAXN
        {
          int N = nNeighUp - GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP;
          int N1 = (N / 4) * 4;
          for (int iUp = 0; iUp < N1; iUp += 4) {
#pragma unroll(4)
            for (int k = 0; k < 4; k++) {
              const float dy = yzdn.x - yzUp[iUp + k];
              const float dz = yzdn.y - yzUp2[iUp + k];
              const float d = dy * dy + dz * dz;
              if (d < bestD) {
                bestD = d;
                bestDn = i;
                bestUp = GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP + iUp + k;
              }
            }
          }
          if (N1 + 1 < N) {
            int iUp = N1 + 1;
            const float dy = yzdn.x - yzUp[iUp];
            const float dz = yzdn.y - yzUp2[iUp];
            const float d = dy * dy + dz * dz;
            if (d < bestD) {
              bestD = d;
              bestDn = i;
              bestUp = GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP + iUp;
            }
          }
          if (N1 + 2 < N) {
            int iUp = N1 + 2;
            const float dy = yzdn.x - yzUp[iUp];
            const float dz = yzdn.y - yzUp2[iUp];
            const float d = dy * dy + dz * dz;
            if (d < bestD) {
              bestD = d;
              bestDn = i;
              bestUp = GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP + iUp;
            }
          }
          if (N1 + 3 < N) {
            int iUp = N1 + 3;
            const float dy = yzdn.x - yzUp[iUp];
            const float dz = yzdn.y - yzUp2[iUp];
            const float d = dy * dy + dz * dz;
            if (d < bestD) {
              bestD = d;
              bestDn = i;
              bestUp = GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP + iUp;
            }
          }
        }
#endif
      }
    }

    if (bestD <= chi2Cut) {
#if GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP > 0 && GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP < GPUCA_MAXN
      linkUp = bestUp >= GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP ? neighUp[bestUp - GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP] : s.mB[bestUp][iThread];
#elif GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP == GPUCA_MAXN
      linkUp = s.mB[bestUp][iThread];
#else
      linkUp = neighUp[bestUp];
#endif
      linkDn = bestDn;
    }

    tracker.mData.mLinkUpData[lHitNumberOffset + ih] = linkUp;
    tracker.mData.mLinkDownData[lHitNumberOffset + ih] = linkDn;
  }
}
