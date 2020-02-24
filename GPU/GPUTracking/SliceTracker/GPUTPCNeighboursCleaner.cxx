// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUTPCNeighboursCleaner.cxx
/// \author Sergey Gorbunov, Ivan Kisel, David Rohr

#include "GPUTPCNeighboursCleaner.h"
#include "GPUTPCTracker.h"
#include "GPUCommonMath.h"

using namespace GPUCA_NAMESPACE::gpu;

template <>
GPUdii() void GPUTPCNeighboursCleaner::Thread<0>(int /*nBlocks*/, int nThreads, int iBlock, int iThread, GPUsharedref() MEM_LOCAL(GPUSharedMemory) & s, processorType& tracker)
{
  unsigned int iRow = iBlock + 2;
  unsigned int nHits = tracker.Row(iRow).NHits();
  unsigned int* linkUpData = (unsigned int*)(tracker.Data().mLinkUpData);

#ifdef GPUCA_GPUCODE  
  const GPUTPCRow& __restrict__ row = tracker.Row(iRow);
#else
  const GPUTPCRow& row = tracker.Row(iRow);
#endif

 const unsigned int &rowOffset = row.mHitNumberOffset;

  for (long iter = 0; iter < 1000; iter++) {
    for (unsigned int ih = iThread; ih < nHits; ih += nThreads) {
      linkUpData[row.mHitNumberOffset + ih] = CALINK_INVAL;
      //linkUpData[rowOffset + ih] = CALINK_INVAL;
    }
  }
}
