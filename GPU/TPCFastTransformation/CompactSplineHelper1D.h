// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  CompactSplineHelper1D.h
/// \brief Definition of CompactSplineHelper1D class

/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#ifndef ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_COMPACTSPLINEHELPER1D_H
#define ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_COMPACTSPLINEHELPER1D_H

#include <cmath>
#include <vector>

#include "GPUCommonDef.h"
#include "Rtypes.h"
#include "TString.h"
#include "CompactSpline1D.h"
#include <functional>

namespace GPUCA_NAMESPACE
{
namespace gpu
{

///
/// The CompactSplineHelper1D class is to initialize CompactSpline* objects
///

class CompactSplineHelper1D
{
 public:
  ///
  /// \brief Helper structure for 1D spline construction
  ///
  struct Point {
    int iKnot; ///< knot index
    double cf0;
    double cz0;
    double cf1;
    double cz1;
  };

  /// _____________  Constructors / destructors __________________________

  /// Default constructor
  CompactSplineHelper1D();

  /// Copy constructor: disabled
  CompactSplineHelper1D(const CompactSplineHelper1D&) CON_DELETE;

  /// Assignment operator: disabled
  CompactSplineHelper1D& operator=(const CompactSplineHelper1D&) CON_DELETE;

  /// Destructor
  ~CompactSplineHelper1D() CON_DEFAULT;

  /// _______________  Main functionality  ________________________

  /// Creates classical spline data for a given input function
  std::unique_ptr<float[]> constructDataClassical1D(const CompactSpline1D& spline, std::function<float(float)> F, float uMin, float uMax);

  /// Creates compact spline data for a given input function
  std::unique_ptr<float[]> constructData1D(const CompactSpline1D& spline, std::function<float(float)> F, float uMin, float uMax, int nAxiliaryPoints);

  /// _______________   Tools for a manual construction of compact splines   ________________________

  int setSpline(const CompactSpline1D& spline, int nAxiliaryPoints);
  int getNdataPoints() const { return mPoints.size(); }

  /// N parameters in the data array per output dimension
  int getNparameters() const { return 2 * mNKnots; }

  void constructData1D(const float inF[/*getNdataPoints()*/], float outSplineData[/*getNparameters()*/]) const;
  void constructDataGradually(int Ndim, const float inF[/*N Data Points x Ndim */], float outSplineData[/*N Spline Parameters*/]) const;

  template <int Ndim>
  void constructDataGradually(const float inF[/*N Data Points x Ndim */], float outSplineData[/*N Spline Parameters*/]) const
  {
    constructDataGradually(Ndim, inF, outSplineData);
  }
  /// _______________  Utilities   ________________________

  int getKnotPoint(int iknot) const { return mKnotPoints[iknot]; }

  ///  Gives error string
  const char* getLastError() const { return mError.Data(); }

 private:
  /// Stores an error message
  int storeError(Int_t code, const char* msg);

  TString mError = ""; ///< error string

  /// helpers for the construction of 1D spline
  int mNKnots = 0;
  std::vector<Point> mPoints;
  std::vector<int> mKnotPoints;
  std::vector<double> mMatrixI;
  std::vector<double> mMatrixFastI;
  std::vector<double> mMatrixFastF;
};

inline int CompactSplineHelper1D::storeError(int code, const char* msg)
{
  mError = msg;
  return code;
}
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
