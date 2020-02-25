// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUTPCCFDecodeZS.cxx
/// \author David Rohr

#include "GPUTPCCFDecodeZS.h"
#include "GPUCommonMath.h"
#include "GPUTPCClusterFinder.h"
#include "DataFormatsTPC/ZeroSuppression.h"

#ifndef __OPENCL__
#include "Headers/RAWDataHeader.h"
#else
namespace o2
{
namespace header
{
struct RAWDataHeader {
  unsigned int words[16];
};
} // namespace header
} // namespace o2

#endif

using namespace GPUCA_NAMESPACE::gpu;
using namespace o2::tpc;

template <>
GPUdii() void GPUTPCCFDecodeZS::Thread<GPUTPCCFDecodeZS::decodeZS>(int nBlocks, int nThreads, int iBlock, int iThread, GPUSharedMemory& smem, processorType& clusterer)
{
  GPUTPCCFDecodeZS::decode(clusterer, smem, nBlocks, nThreads, iBlock, iThread);
}

GPUd() void GPUTPCCFDecodeZS::decode(GPUTPCClusterFinder& clusterer, GPUSharedMemory& s, int nBlocks, int nThreads, int iBlock, int iThread)
{
  const unsigned int slice = clusterer.mISlice;
  const unsigned int endpoint = iBlock;
  const GPUTrackingInOutZS::GPUTrackingInOutZSSlice& zs = clusterer.GetConstantMem()->ioPtrs.tpcZS->slice[slice];

  //unsigned int s_RowClusterOffset[o2::tpc::TPCZSHDR::TPC_MAX_ZS_ROW_IN_ENDPOINT];

  if (zs.count[endpoint] == 0) {
    return;
  }
  deprecated::PackedDigit* digits = clusterer.mPdigits;
  const size_t nDigits = clusterer.mPmemory->nDigitsOffset[endpoint];
  unsigned int rowOffsetCounter = 0;

  const unsigned int s_decodeBits = TPCZSHDR::TPC_ZS_NBITS_V1;

  unsigned int tmpOutput = 0;

  const unsigned int* pageSrc = (const unsigned int*)(((const unsigned char*)zs.zsPtr[endpoint][0]));

  GPUbarrier();
  CA_SHARED_CACHE_REF(&s.ZSPage[0], pageSrc, TPCZSHDR::TPC_ZS_PAGE_SIZE, unsigned int, pageCache);
  GPUbarrier();

  if (iThread != 0)
    return;

  for (int iter = 0; iter < 100000; iter++) {
    const unsigned char* __restrict__ page = (const unsigned char*)pageCache;
    const unsigned char* pagePtr = page + sizeof(o2::header::RAWDataHeader);
    const TPCZSHDR* hdr = reinterpret_cast<const TPCZSHDR*>(pagePtr);
    pagePtr += sizeof(*hdr);
    const int nTimeBins = hdr->nTimeBins;
    //for (int l = 0; l < nTimeBins; l++) {
         for (int l = 0; l < 1; l++) {
 
      pagePtr += (pagePtr - page) & 1; //Ensure 16 bit alignment
      const TPCZSTBHDR* tbHdr = reinterpret_cast<const TPCZSTBHDR*>(pagePtr);
      const unsigned int mask = ((const unsigned int)tbHdr->rowMask) & 0x7FFF;
      if (mask == 0) {
        pagePtr += 2;
        continue;
      }
      const int nRowsUsed = CAMath::Popcount(mask);
      pagePtr += 2 * nRowsUsed;
      const auto sg1 = tbHdr->rowAddr1();
      //GPUbarrier();
      { //if (iThread == 0) {
        {
          //s_RowClusterOffset[0] = rowOffsetCounter;
          const unsigned char* __restrict__ rowData = pagePtr;
          unsigned int tmp = *rowData;
          rowOffsetCounter += rowData[tmp << 1]; // Sum up number of ADC samples per row to compute offset in target buffer
        }
        for (int n = 1; n < nRowsUsed; n++) {
          //s_RowClusterOffset[n] = rowOffsetCounter;
          const unsigned char* __restrict__ rowData = (page + sg1[n - 1]);
          unsigned int tmp = *rowData;
          rowOffsetCounter += rowData[tmp << 1]; // Sum up number of ADC samples per row to compute offset in target buffer
        }
      }
      //GPUbarrier();
      tmpOutput += rowOffsetCounter;
      if (nRowsUsed > 1) {
        //pagePtr = page + tbHdr->rowAddr1()[nRowsUsed - 2];
        pagePtr = page + sg1[nRowsUsed - 2];
      }
      pagePtr += 2 * *pagePtr;                          // Go to entry for last sequence length
      pagePtr += 1 + (*pagePtr * s_decodeBits + 7) / 8; // Go to beginning of next time bin
    }
  }

  if (iThread == 0) {
    digits[nDigits].time = tmpOutput;
  }
}
