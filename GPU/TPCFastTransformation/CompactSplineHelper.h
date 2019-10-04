// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  CompactSplineHelper.h
/// \brief Definition of CompactSplineHelper class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#ifndef ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_COMPACTSPLINEHELPER_H
#define ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_COMPACTSPLINEHELPER_H

#include <cmath>

#include "GPUCommonDef.h"
#include "Rtypes.h"
#include "TString.h"
#include "CompactSplineIrregular1D.h"
#include "CompactSplineIrregular2D3D.h"
#include <functional>

namespace GPUCA_NAMESPACE
{
namespace gpu
{

///
/// The CompactSplineHelper class is to initialize CompactSpline* objects
///

class CompactSplineHelper
{
 public:
  /// _____________  Constructors / destructors __________________________

  /// Default constructor
  CompactSplineHelper();

  /// Copy constructor: disabled
  CompactSplineHelper(const CompactSplineHelper&) CON_DELETE;

  /// Assignment operator: disabled
  CompactSplineHelper& operator=(const CompactSplineHelper&) CON_DELETE;

  /// Destructor
  ~CompactSplineHelper() CON_DEFAULT;

  /// _______________  Main functionality  ________________________

  /// Creates spline data for a given input function
  std::unique_ptr<float[]> createClassical(const CompactSplineIrregular1D& spline, std::function<float(float)> F);

  std::unique_ptr<float[]> create(const CompactSplineIrregular1D& spline, const double inputU[], const double inputF[], int inputN);
  std::unique_ptr<float[]> create(const CompactSplineIrregular1D& spline, std::function<double(double)> F, int nAxiliaryPoints = 2);

  std::unique_ptr<float[]> create(const CompactSplineIrregular2D3D& spline, std::function<void(float, float, float&, float&, float&)> F, int nAxiliaryPoints = 2);

  /// _______________  Utilities   ________________________

  ///  Gives error string
  const char* getLastError() const { return mError.Data(); }

 private:
  /// Stores an error message
  int storeError(Int_t code, const char* msg);

  TString mError; ///< error string
};

inline int CompactSplineHelper::storeError(int code, const char* msg)
{
  mError = msg;
  return code;
}
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
