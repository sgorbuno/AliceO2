// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  SplineHelper1D.h
/// \brief Definition of SplineHelper1D class

/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#ifndef ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SPLINEHELPER1D_H
#define ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SPLINEHELPER1D_H

#include <cmath>
#include <vector>

#include "GPUCommonDef.h"
#include "Rtypes.h"
#include "TString.h"
#include "Spline1D.h"
#include <functional>

namespace GPUCA_NAMESPACE
{
namespace gpu
{

///
/// The SplineHelper1D class is to initialize Spline* objects
///

class SplineHelper1D
{
 public:
  ///
  /// \brief Helper structure for 1D spline construction
  ///
  struct MeasurementPoint {
    double u; ///< u coordinate 
    double cS0; ///< a coefficient for s0
    double cZ0; ///< a coefficient for s'0
    double cS1; ///< a coefficient for s1
    double cZ1; ///< a coefficient for s'1
    int iKnot; ///< index of the left knot of the segment
    bool isKnot; ///< is the point placed at a knot
  };

  /// _____________  Constructors / destructors __________________________

  /// Default constructor
  SplineHelper1D();

  /// Copy constructor: disabled
  SplineHelper1D(const SplineHelper1D&) CON_DELETE;

  /// Assignment operator: disabled
  SplineHelper1D& operator=(const SplineHelper1D&) CON_DELETE;

  /// Destructor
  ~SplineHelper1D() CON_DEFAULT;

  /// _______________  Main functionality  ________________________

  int setSpline(const Spline1D& spline, int nAxiliaryPoints);

  /// Create classical spline parameters for a given input function
  std::unique_ptr<float[]> constructParametersClassical(int Ndim, std::function<void(float, float[])> F, float uMin, float uMax);

  /// Create compact spline parameters for a given input function
  std::unique_ptr<float[]> constructParameters(int Ndim, std::function<void(float, float[])> F, float uMin, float uMax);

  /// Create compact spline parameters gradually
  std::unique_ptr<float[]> constructParametersGradually(int Ndim, std::function<void(float, float[])> F, float uMin, float uMax);

  /// _______________   Interface for a manual construction of compact splines   ________________________

  int getNumberOfMeasurements() const { return mMeasurementPoints.size(); }

  void constructParameters(int Ndim, const float F[/*getNumberOfMeasurements() x Ndim*/], float parameters[/*mSpline.getNumberOfParameters(Ndim)*/]) const;

  void constructParametersGradually(int Ndim, const float F[/*getNumberOfMeasurements() x Ndim */], float parameters[/*mSpline.getNumberOfParameters(Ndim)*/]) const;

  void copySfromMeasurements(int Ndim, const float F[/*getNumberOfMeasurements() x Ndim*/], float parameters[/*mSpline.getNumberOfParameters(Ndim)*/]) const;

  void constructDerivatives(int Ndim, const float F[/*getNumberOfMeasurements() x Ndim*/], float parameters[/*mSpline.getNumberOfParameters(Ndim)*/]) const;

  /// _______________  Utilities   ________________________

  const Spline1D& getSpline() const { return mSpline; }

  int getKnotMeasurement(int iknot) const { return mKnotMeasurements[iknot]; }

  const MeasurementPoint& getMeasurementPoint(int ip) const { return mMeasurementPoints[ip]; }

  ///  Gives error string
  const char* getLastError() const { return mError.Data(); }

 private:
  /// Stores an error message
  int storeError(Int_t code, const char* msg);

  TString mError = ""; ///< error string

  /// helpers for the construction of 1D spline

  Spline1D mSpline; ///< copy of the spline
  std::vector<MeasurementPoint> mMeasurementPoints; ///< measurement points
  std::vector<int> mKnotMeasurements; ///< which measurement points are at knots
  std::vector<double> mLSMmatrixFull; ///< a matrix to convert the measurements into the spline parameters with the LSM method
  std::vector<double> mLSMmatrixSderivatives;
  std::vector<double> mLSMmatrixSvalues;
};

inline int SplineHelper1D::storeError(int code, const char* msg)
{
  mError = msg;
  return code;
}
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
