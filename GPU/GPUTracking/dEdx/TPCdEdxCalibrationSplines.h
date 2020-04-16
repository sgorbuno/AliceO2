// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  dEdxCalibrationSplines.h
/// \brief Definition of dEdxCalibrationSplines class
///
/// \author  Matthias Kleiner <matthias.kleiner@cern.ch>

#ifndef TPCdEdxCalibrationSplines_H
#define TPCdEdxCalibrationSplines_H

#include "FlatObject.h"
#include "Spline2D.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{

///
/// The dEdxCalibrationSplines class represents the calibration of the dEdx of mostly geometrical effects
///

class TPCdEdxCalibrationSplines : public FlatObject
{
 public:
  typedef Spline2D<float, 1, 1> SplineType;

  /// _____________  Constructors / destructors __________________________

  /// Default constructor: creates an empty uninitialized object
  TPCdEdxCalibrationSplines();

  /// Copy constructor: disabled to avoid ambiguity.
  TPCdEdxCalibrationSplines(const TPCdEdxCalibrationSplines&) CON_DELETE;

  /// Assignment operator: disabled to avoid ambiguity.
  TPCdEdxCalibrationSplines& operator=(const TPCdEdxCalibrationSplines&) CON_DELETE;

  /// Destructor
  ~TPCdEdxCalibrationSplines();

  /// Gives pointer to a spline
  GPUd() const SplineType& getSpline(int chargeType, int region) const;

  /// sets the splines from an input file
  void setSplinesFromFile(TFile& inpf);

  /// returns the number of splines stored in the calibration object
  GPUd() unsigned int getFSplines() const
  {
    return mFSplines;
  };

  GPUd() float interpolateqMax(const int splineInd, const float angleZ, const float z) const
  {
    return mCalibSplinesqMax[splineInd].interpolate(angleZ, z);
  };

  GPUd() float interpolateqTot(const int splineInd, const float angleZ, const float z) const
  {
    return mCalibSplinesqTot[splineInd].interpolate(angleZ, z);
  };

  GPUd() SplineType& getSplineqMax(const int splineInd)
  {
    return mCalibSplinesqMax[splineInd];
  };

  GPUd() SplineType& getSplineqTot(const int splineInd)
  {
    return mCalibSplinesqTot[splineInd];
  };

    /// _______________  IO   ________________________
#if !defined(GPUCA_ALIGPUCODE) && !defined(GPUCA_STANDALONE)
  /// read a class object from the file
  static TPCdEdxCalibrationSplines* readFromFile(TFile& inpf, const char* name);
#endif

 private:
  constexpr static unsigned int mFSplines = 10; ///< number of splines stored for each type
  SplineType mCalibSplinesqMax[mFSplines];      ///< spline objects storage for the splines for qMax
  SplineType mCalibSplinesqTot[mFSplines];      ///< spline objects storage for the splines for qTot

  ClassDefNV(TPCdEdxCalibrationSplines, 1);
};

#if !defined(GPUCA_ALIGPUCODE) && !defined(GPUCA_STANDALONE)

TPCdEdxCalibrationSplines* TPCdEdxCalibrationSplines::readFromFile(
  TFile& inpf, const char* name)
{
  /// read a class object from the file
  return FlatObject::readFromFile<TPCdEdxCalibrationSplines>(inpf, name);
}

void TPCdEdxCalibrationSplines::setSplinesFromFile(TFile& inpf)
{
  for (int ireg = 0; ireg <= 9; ++ireg) {
    o2::gpu::Spline2D<float, 1>* splineTmpqMax = o2::gpu::Spline2D<float, 1>::readFromFile(inpf, Form("spline_qMax_region%d", ireg));
    o2::gpu::Spline2D<float, 1>* splineTmpqTot = o2::gpu::Spline2D<float, 1>::readFromFile(inpf, Form("spline_qTot_region%d", ireg));
    mCalibSplinesqMax[ireg] = *splineTmpqMax;
    mCalibSplinesqTot[ireg] = *splineTmpqTot;
  }
}

#endif

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
