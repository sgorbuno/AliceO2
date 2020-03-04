
#include "GPUTPCNeighboursFinder.h"
using namespace GPUCA_NAMESPACE::gpu;

template <>
GPUdii() void GPUTPCNeighboursFinder::Thread<0>(int /*nBlocks*/, int nThreads, int iBlock, int iThread, GPUsharedref() MEM_LOCAL(GPUSharedMemory) & GPUrestrict() ss, processorType& GPUrestrict() trackerX)
{
#ifdef GPUCA_GPUCODE

#define doCompileCrash 1

  /*  Problem 1:  --- Compilation crash --- 

   The code does not not compile when doCompileCrash == 1
   and compiles otherwise.
  
   Compilation error:

   [ 98%] Built target ca
    error: local memory limit exceeded (69916) in _Z27krnl_GPUTPCNeighboursFinderi
    Generating AMD GCN kernel failed in llc for target: gfx906
    Error: hc-kernel-assemble[166]: failed with status -1
    clang-10: error: HC assembler command failed with exit code 255 (use -v to see invocation)
    GPU/GPUTracking/Base/hip/CMakeFiles/GPUTrackingHIP.dir/build.make:62: recipe for target 'GPU/GPUTracking/Base/hip/CMakeFiles/GPUTrackingHIP.dir/GPUReconstructionHIP.hip.cxx.o' failed
    make[2]: *** [GPU/GPUTracking/Base/hip/CMakeFiles/GPUTrackingHIP.dir/GPUReconstructionHIP.hip.cxx.o] Error 255
    CMakeFiles/Makefile2:290: recipe for target 'GPU/GPUTracking/Base/hip/CMakeFiles/GPUTrackingHIP.dir/all' failed
    make[1]: *** [GPU/GPUTracking/Base/hip/CMakeFiles/GPUTrackingHIP.dir/all] Error 2
    Makefile:129: recipe for target 'all' failed
    make: *** [all] Error 2

  */

  /*  Problem 2: 
    ss.mRow.Grid().GetBin(..) is not optimised away, 
    though Grid() and GetBin() are declared const and the output values 
    binYmin, binZmin are not used anywhere.
  */

#if doCompileCrash == 1
// do nothing
#else
  ss.mB[0][0] = 0; // write something to the shared memory
#endif

  int binYmin, binZmin;
  reinterpret_cast<const GPUTPCRow&>(ss.mRow).Grid().GetBin(0.f, 0.f, &binYmin, &binZmin);
#endif
}
