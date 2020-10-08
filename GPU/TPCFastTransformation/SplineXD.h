// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  SplineXD.h
/// \brief Definition of SplineXD class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#ifndef ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SPLINEXD_H
#define ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SPLINEXD_H

#include "GPUCommonDef.h"
#include "FlatObject.h"
#if !defined(GPUCA_GPUCODE)
#include <functional>
#endif

#include <limits>

class TFile;

namespace GPUCA_NAMESPACE
{
namespace gpu
{

class blabla
{
  ClassDefNV(blabla, 1);
};

template <typename DataT, int nYdimT = std::numeric_limits<int>::min()>
class SplineXD;

template <typename DataT>
class SplineXD<DataT, std::numeric_limits<int>::min()>
{
 public:
  SplineXD(int nYdim) : mYdim(nYdim), mX(0.) {}

  int mYdim;
  DataT mX;
  ClassDefNV(SplineXD, 0);
};

template <typename DataT, int nYdimT>
class SplineXD //: public SplineXD<DataT>
{
 public:
  //typedef SplineXD<DataT> TBase;

  /// Default constructor
  SplineXD() : mYdim(nYdimT), mX(1.) {}

  /// Destructor
  ~SplineXD() CON_DEFAULT;

  int mYdim;
  DataT mX;

  //  ClassDefNV(SplineXD, 0);
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
