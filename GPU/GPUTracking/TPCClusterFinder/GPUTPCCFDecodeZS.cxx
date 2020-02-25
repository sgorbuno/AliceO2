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

  if (zs.count[endpoint] == 0) {
    // return;
  }
  deprecated::PackedDigit* digits = clusterer.mPdigits;
  //const size_t nDigits = clusterer.mPmemory->nDigitsOffset[endpoint];
  /*
  const unsigned int* pageSrc = (const unsigned int*)(((const unsigned char*)zs.zsPtr[endpoint][0]));   

   unsigned char* page = ( unsigned char*)pageSrc;
   unsigned char* pagePtr = page + sizeof(o2::header::RAWDataHeader) + sizeof(TPCZSHDR);
  pagePtr += (pagePtr - page) & 1; //Ensure 16 bit alignment
   TPCZSTBHDR* tbHdr = reinterpret_cast< TPCZSTBHDR*>(pagePtr);
  //const int nRowsUsed = CAMath::Popcount(tbHdr->rowMask & 0x7FFF);
*/
GPUshared()  unsigned char pp[1000];
unsigned char* arr = reinterpret_cast< unsigned char*> (pp);

  if (iThread != 0)
    return;
 

for( int i=0; i<9; i++) arr[i] = 10*i;
pp[0]=9;

GPUbarrier();

const int nRowsUsed = pp[0];

  unsigned int tmpOutput = 0;

  for (int iter = 0; iter < 1000000; iter++) {
    for (int n = 1; n < nRowsUsed; n++) {
      const unsigned char* rowData = (pp + arr[n - 1]);
      tmpOutput += rowData[2 * *rowData];
    }
  }

  if (iThread == 0) {
    digits[0].time = tmpOutput;
  }
}
