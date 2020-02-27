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
using namespace GPUCA_NAMESPACE::gpu;

#if defined(__HIPCC__) && defined(GPUCA_GPUCODE_DEVICE)
#define HIPGPUsharedref() __attribute__((address_space(3)))
#define HIPGPUglobalref() __attribute__((address_space(1)))
#define HIPGPUconstantref() __attribute__((address_space(4)))
#else
#define HIPGPUsharedref()
#define HIPGPUglobalref()
#define HIPGPUconstantref()
#endif

template <>
GPUdii() void GPUTPCNeighboursFinder::Thread<0>(int /*nBlocks*/, int nThreads, int iBlock, int iThread, GPUsharedref() MEM_LOCAL(GPUSharedMemory) & GPUrestrict() ss, processorType& GPUrestrict() trackerX)
{
  //* find neighbours

#ifdef GPUCA_GPUCODE
  for (unsigned int i = iThread; i < sizeof(MEM_PLAIN(GPUTPCRow)) / sizeof(int); i += nThreads) {
    ((int*)&ss.mRow)[i] = ((int*)&trackerX.SliceDataRows()[iBlock])[i];
    if (iBlock >= 2 && iBlock < GPUCA_ROW_COUNT - 2) {
      ((int*)&ss.mRowUp)[i] = ((int*)&trackerX.SliceDataRows()[iBlock + 2])[i];
      ((int*)&ss.mRowDown)[i] = ((int*)&trackerX.SliceDataRows()[iBlock - 2])[i];
    }
  }
  GPUbarrier();
#endif

  HIPGPUsharedref() MEM_LOCAL(GPUSharedMemory) & GPUrestrict() s = (HIPGPUsharedref() MEM_LOCAL(GPUSharedMemory)&)ss;
  HIPGPUconstantref() processorType& GPUrestrict() tracker = (HIPGPUconstantref() processorType&)trackerX;

  if (iThread == 0) {
    s.mIRow = iBlock;
    if (s.mIRow < GPUCA_ROW_COUNT) {
#ifdef GPUCA_GPUCODE
      HIPGPUsharedref() const MEM_LOCAL(GPUTPCRow) & GPUrestrict() row = s.mRow;
#else
      HIPGPUglobalref() const MEM_GLOBAL(GPUTPCRow) & GPUrestrict() row = tracker.mData.mRows[s.mIRow];
#endif
      s.mNHits = row.mNHits;

      if ((s.mIRow >= 2) && (s.mIRow <= GPUCA_ROW_COUNT - 3)) {
        s.mIRowUp = s.mIRow + 2;
        s.mIRowDn = s.mIRow - 2;

        // references to the rows above and below

#ifdef GPUCA_GPUCODE
        HIPGPUsharedref() const MEM_LOCAL(GPUTPCRow) & GPUrestrict() rowUp = s.mRowUp;
        HIPGPUsharedref() const MEM_LOCAL(GPUTPCRow) & GPUrestrict() rowDn = s.mRowDown;
#else
        HIPGPUglobalref() const MEM_GLOBAL(GPUTPCRow) & GPUrestrict() rowUp = tracker.mData.mRows[s.mIRowUp];
        HIPGPUglobalref() const MEM_GLOBAL(GPUTPCRow) & GPUrestrict() rowDn = tracker.mData.mRows[s.mIRowDn];
#endif
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

  if ((s.mIRow <= 1) || (s.mIRow >= GPUCA_ROW_COUNT - 2)) {
    return;
  }

  const float chi2Cut = 3.f * 3.f * 4 * (s.mUpDx * s.mUpDx + s.mDnDx * s.mDnDx);
  const int chi2Tmp = (int)chi2Cut;
// float chi2Cut = 3.*3.*(s.mUpDx*s.mUpDx + s.mDnDx*s.mDnDx ); //SG
#ifdef GPUCA_GPUCODE
  HIPGPUsharedref() const MEM_LOCAL(GPUTPCRow) & GPUrestrict() row = (HIPGPUsharedref() const MEM_LOCAL(GPUTPCRow)&)s.mRow;
  HIPGPUsharedref() const MEM_LOCAL(GPUTPCRow) & GPUrestrict() rowUp = (HIPGPUsharedref() const MEM_LOCAL(GPUTPCRow)&)s.mRowUp;
  HIPGPUsharedref() const MEM_LOCAL(GPUTPCRow) & GPUrestrict() rowDn = (HIPGPUsharedref() const MEM_LOCAL(GPUTPCRow)&)s.mRowDown;
#else
  HIPGPUglobalref() const MEM_GLOBAL(GPUTPCRow) & GPUrestrict() row = (HIPGPUglobalref() const MEM_GLOBAL(GPUTPCRow)&)tracker.mData.mRows[s.mIRow];
  HIPGPUglobalref() const MEM_GLOBAL(GPUTPCRow) & GPUrestrict() rowUp = (HIPGPUglobalref() const MEM_GLOBAL(GPUTPCRow)&)tracker.mData.mRows[s.mIRowUp];
  HIPGPUglobalref() const MEM_GLOBAL(GPUTPCRow) & GPUrestrict() rowDn = (HIPGPUglobalref() const MEM_GLOBAL(GPUTPCRow)&)tracker.mData.mRows[s.mIRowDn];
#endif
  const float y0 = row.mGrid.mYMin;
  const float z0 = row.mGrid.mZMin;
  const float stepY = row.mHstepY;
  const float stepZ = row.mHstepZ;

  const long int lHitNumberOffset = row.mHitNumberOffset;
  const long int lHitNumberOffsetDn = rowDn.mHitNumberOffset;  
  const int lFirstHitInBinOffsetDn = rowDn.mFirstHitInBinOffset;
  HIPGPUglobalref() const calink* GPUrestrict() lFirstHitInBin = (HIPGPUglobalref() const calink*)tracker.mData.mFirstHitInBin;
  HIPGPUglobalref() const cahit2* GPUrestrict() pHitData = (HIPGPUglobalref() const cahit2*)tracker.mData.mHitData;

  const float y0Dn = rowDn.mGrid.mYMin;
  const float z0Dn = rowDn.mGrid.mZMin;
  const float stepYDn = rowDn.mHstepY;
  const float stepZDn = rowDn.mHstepZ;

  for (int ih = iThread; ih < s.mNHits; ih += nThreads) {

    int linkDn = -1;

    if (rowDn.mNHits > 0) {
      HIPGPUglobalref() const cahit2& hitData = pHitData[lHitNumberOffset + ih];
      const float y = y0 + (hitData.x) * stepY;
      const float z = z0 + (hitData.y) * stepZ;

      float minZ, maxZ, minY, maxY;
      int binYmin, binYmax, binZmin, binZmax;
      int nY;

      const float kAngularMultiplier = tracker.mConstantMem->param.rec.SearchWindowDZDR;
      const float kAreaSize = tracker.mConstantMem->param.rec.NeighboursSearchArea;


#if GPUCA_NEIGHBOURS_FINDER_MAX_NNEIGHUP > 0
          s.mB[0][iThread] = (calink)0;        
#endif

      if (1) {
        {
          const float yy = y * s.mDnTx;
          const float zz = kAngularMultiplier != 0.f ? z : (z * s.mDnTx);
          const float dy = kAreaSize;
          const float dz = kAngularMultiplier != 0.f ? (-s.mDnDx * kAngularMultiplier) : kAreaSize;
          minZ = zz - dz;
          maxZ = zz + dz;
          minY = yy - dy;
          maxY = yy + dy;
          reinterpret_cast<const GPUTPCRow&>(rowDn).Grid().GetBin(minY, minZ, &binYmin, &binZmin);
          reinterpret_cast<const GPUTPCRow&>(rowDn).Grid().GetBin(maxY, maxZ, &binYmax, &binZmax);
          nY = reinterpret_cast<const GPUTPCRow&>(rowDn).Grid().Ny();
        }
        int bestDn = -1, bestUp = -1;
        float bestD = 1.e10f;

       // int nNeighDn = 0;
        for (int k1 = binZmin; k1 <= binZmax; k1++) {
          int iMin = lFirstHitInBin[lFirstHitInBinOffsetDn + k1 * nY + binYmin];
          int iMax = lFirstHitInBin[lFirstHitInBinOffsetDn + k1 * nY + binYmax + 1];
          for (int i = iMin; i < iMax; i++) {
            HIPGPUglobalref() const cahit2& hitDataDn = pHitData[lHitNumberOffsetDn + i];
            float d = hitDataDn.x + hitDataDn.y;
            /*
            GPUTPCHit h;
            h.mY = y0Dn + (hitDataDn.x) * stepYDn;
            h.mZ = z0Dn + (hitDataDn.y) * stepZDn;
            */
            /*
            if (h.mY < minY || h.mY > maxY || h.mZ < minZ || h.mZ > maxZ)
              continue;

            nNeighDn++;
*/
            if (d < bestD) {
              bestD = d;
              bestDn = i;              
            }
          }
        }

        if (1) {
          linkDn = bestDn;
        }
      }
    }

    tracker.mData.mLinkDownData[lHitNumberOffset + ih] = chi2Tmp + linkDn;
  }
}
