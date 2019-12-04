// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  CompactSplineHelper2D.h
/// \brief Definition of CompactSplineHelper2D class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#ifndef ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_COMPACTSPLINEHELPER2D_H
#define ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_COMPACTSPLINEHELPER2D_H

#include <cmath>
#include <vector>

#include "GPUCommonDef.h"
#include "Rtypes.h"
#include "TString.h"
#include "CompactSpline1D.h"
#include "CompactSpline2D.h"
#include "CompactSplineHelper1D.h"
#include <functional>

namespace GPUCA_NAMESPACE
{
namespace gpu
{

///
/// The CompactSplineHelper2D class is to initialize CompactSpline* objects
///

class CompactSplineHelper2D
{
 public:
  /// _____________  Constructors / destructors __________________________

  /// Default constructor
  CompactSplineHelper2D();

  /// Copy constructor: disabled
  CompactSplineHelper2D(const CompactSplineHelper2D&) CON_DELETE;

  /// Assignment operator: disabled
  CompactSplineHelper2D& operator=(const CompactSplineHelper2D&) CON_DELETE;

  /// Destructor
  ~CompactSplineHelper2D() CON_DEFAULT;

  /// _______________  Main functionality  ________________________

  /// Creates compact spline data for a given input function
  //std::unique_ptr<float[]> constructSpline2D(const CompactSpline2DD& spline, std::function<float(float)> F, float uMin, float uMax, int nAxiliaryPoints);

  /// Tools for a manual construction of compact splines
  int setSpline(const CompactSpline2D& spline, int nAxiliaryPointsU, int nAxiliaryPointsV);
  int getNdataPointsU() const { return mHelperU.getNdataPoints(); }
  int getNdataPointsV() const { return mHelperV.getNdataPoints(); }
  int getNdataPoints() const { return getNdataPointsU() * getNdataPointsV(); }
  int getNparameters() const { return mSpline.getDataSizeInElements<1>(); }

  void constructSpline(const float inF[/*getNdataPoints()*/], float outSplineData[/*getNparameters()*/]) const;

  /// _______________  Utilities   ________________________

  ///  Gives error string
  const char* getLastError() const { return mError.Data(); }

 private:
  /// Stores an error message
  int storeError(Int_t code, const char* msg);

  TString mError = ""; ///< error string

  CompactSpline2D mSpline;
  CompactSplineHelper1D mHelperU;
  CompactSplineHelper1D mHelperV;
};

inline int CompactSplineHelper2D::storeError(int code, const char* msg)
{
  mError = msg;
  return code;
}
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
