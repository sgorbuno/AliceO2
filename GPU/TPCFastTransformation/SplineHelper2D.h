// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  SplineHelper2D.h
/// \brief Definition of SplineHelper2D class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#ifndef ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SPLINEHELPER2D_H
#define ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SPLINEHELPER2D_H

#include <cmath>
#include <vector>

#include "GPUCommonDef.h"
#include "GPUCommonRtypes.h"
#include "Spline1D.h"
#include "Spline2D.h"
#include "SplineHelper1D.h"
#include <functional>
#include <string>

namespace GPUCA_NAMESPACE
{
namespace gpu
{

///
/// The SplineHelper2D class is to initialize Spline* objects
///
template <typename DataT>
class SplineHelper2D
{
 public:
  /// _____________  Constructors / destructors __________________________

  /// Default constructor
  SplineHelper2D();

  /// Copy constructor: disabled
  SplineHelper2D(const SplineHelper2D&) CON_DELETE;

  /// Assignment operator: disabled
  SplineHelper2D& operator=(const SplineHelper2D&) CON_DELETE;

  /// Destructor
  ~SplineHelper2D() CON_DEFAULT;

  /// _______________  Main functionality  ________________________

  /// Create best-fit spline parameters for a given input function F
  void approximateFunction(
    Spline2D<DataT>& spline,
    double x1Min, double x1Max, double x2Min, double x2Max,
    std::function<void(double x1, double x2, double f[/*spline.getYdimensions()*/])> F,
    int nAuxiliaryDataPointsU1 = 4, int nAuxiliaryDataPointsU2 = 4);

  /// _______________   Interface for a step-wise construction of the best-fit spline   ________________________

  /// precompute everything needed for the construction
  int setSpline(const Spline2D<DataT>& spline, int nAuxiliaryPointsU1, int nAuxiliaryPointsU2);

  /// approximate std::function, output in Fparameters
  void approximateFunction(
    DataT* Fparameters, double x1Min, double x1Max, double x2Min, double x2Max,
    std::function<void(double x1, double x2, double f[/*mFdimensions*/])> F) const;

  /// approximate std::function, output in Fparameters. F calculates values for a batch of points.
  void approximateFunctionBatch(
    DataT* Fparameters, double x1Min, double x1Max, double x2Min, double x2Max,
    std::function<void(const std::vector<double>& x1, const std::vector<double>& x2, std::vector<double> f[/*mFdimensions*/])> F,
    unsigned int batchsize) const;

  /// approximate a function given as an array of values at data points
  void approximateFunction(
    DataT* Fparameters, const double DataPointF[/*getNumberOfDataPoints() x nFdim*/]) const;

  int getNumberOfDataPointsU1() const { return mHelperU1.getNumberOfDataPoints(); }

  int getNumberOfDataPointsU2() const { return mHelperU2.getNumberOfDataPoints(); }

  int getNumberOfDataPoints() const { return getNumberOfDataPointsU1() * getNumberOfDataPointsU2(); }

  const SplineHelper1D<DataT>& getHelperU1() const { return mHelperU1; }
  const SplineHelper1D<DataT>& getHelperU2() const { return mHelperU2; }

  /// _______________  Utilities   ________________________

  ///  Gives error string
  const char* getLastError() const { return mError.c_str(); }

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE) // code invisible on GPU and in the standalone compilation
  /// Test the Spline2D class functionality
  static int test(const bool draw = 0, const bool drawDataPoints = 1);
#endif

 private:
  /// Stores an error message
  int storeError(int code, const char* msg);

  std::string mError = ""; ///< error string
  int mFdimensions;        ///< n of F dimensions
  SplineHelper1D<DataT> mHelperU1;
  SplineHelper1D<DataT> mHelperU2;

#ifndef GPUCA_ALIROOT_LIB
  ClassDefNV(SplineHelper2D, 0);
#endif
};

template <typename DataT>
void SplineHelper2D<DataT>::approximateFunction(
  Spline2D<DataT>& spline,
  double x1Min, double x1Max, double x2Min, double x2Max,
  std::function<void(double x1, double x2, double f[/*spline.getYdimensions()*/])> F,
  int nAuxiliaryDataPointsU1, int nAuxiliaryDataPointsU2)
{
  /// Create best-fit spline parameters for a given input function F
  setSpline(spline, nAuxiliaryDataPointsU1, nAuxiliaryDataPointsU2);
  approximateFunction(spline.getParameters(), x1Min, x1Max, x2Min, x2Max, F);
  spline.setXrange(x1Min, x1Max, x2Min, x2Max);
}

template <typename DataT>
int SplineHelper2D<DataT>::setSpline(
  const Spline2D<DataT>& spline, int nAuxiliaryPointsU, int nAuxiliaryPointsV)
{
  // Prepare creation of 2D irregular spline
  // The should be at least one (better, two) Auxiliary measurements on each segnment between two knots and at least 2*nKnots measurements in total
  // Returns 0 when the spline can not be constructed with the given nAuxiliaryPoints

  int ret = 0;
  mFdimensions = spline.getYdimensions();
  if (mHelperU1.setSpline(spline.getGridX1(), mFdimensions, nAuxiliaryPointsU) != 0) {
    ret = storeError(-2, "SplineHelper2D::setSpline2D: error by setting U axis");
  }
  if (mHelperU2.setSpline(spline.getGridX2(), mFdimensions, nAuxiliaryPointsV) != 0) {
    ret = storeError(-3, "SplineHelper2D::setSpline2D: error by setting V axis");
  }
  return ret;
}

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
