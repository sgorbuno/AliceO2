// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  IrregularSpline2D3DCalibrator.h
/// \brief Definition of IrregularSpline2D3DCalibrator class
///
/// \author  Oscar Lange
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#ifndef ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SEMIREGULARSPLINE2D3DCALIBRATOR_H
#define ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SEMIREGULARSPLINE2D3DCALIBRATOR_H

#include "GPUCommonDef.h"
#include "IrregularSpline2D3D.h"
#include <memory>
#include <list>

namespace GPUCA_NAMESPACE
{
namespace gpu
{

class SemiregularSpline2D3DCalibrator
{
 public:

  /// _____________  Constructors / destructors __________________________

  /// Default constructor.
  SemiregularSpline2D3DCalibrator();

  /// Destructor
  ~SemiregularSpline2D3DCalibrator() CON_DEFAULT;

  /// set size of the raster grid
  void setRasterSize(int nKnotsU, int nKnotsV);

  /// set maximal size of the spline grid
  void setMaxNKnots(int nKnotsU, int nKnotsV);

  /// set maximal tolerated deviation between the spline and the input function
  void setMaximalDeviation(float maxDeviation)
  {
    mMaxDeviation = maxDeviation;
  }

  /// function constructs a calibrated Spline
  /// \param spline_uv - calibrated spline
  /// \param         F - input function

  std::unique_ptr<float[]> calibrateSpline(SemiregularSpline2D3D& spline_uv, void (*F)(float, float, float&, float&, float&));

  // some getters and step-by-step calibration methods. Only for debugging.

  const SemiregularSpline2D3D& getRaster() const
  {
    return mRaster;
  }
  const float* getRasterData() const
  {
    return mRasterData.data();
  }

  const SemiregularSpline2D3D& getSpline() const
  {
    return mSpline;
  }

  const float* getSplineData() const
  {
    return mSplineData.data();
  }

  void startCalibration(void (*F)(float, float, float&, float&, float&));
  bool doCalibrationStep();

  static constexpr int MinNKnots = 5;
  
 private:
  /// Methods

  void createCurrentSpline();
  void createTrySpline();
  void createSpline(SemiregularSpline2D3D& sp, std::vector<float>& data);

  void getCost(const SemiregularSpline2D3D& spline, const std::vector<float>& data,
	       double &cost, double &maxDeviation) const;

  /// Class members

  int mMaxNKnots[2] = { 5, 5 }; ///< max N knots, U / V axis
  
  int mNKnotsV; ///< number of knots V axis
  std::vector<int> mNKnotsU; ///< vector of knots for U axis per V row

  SemiregularSpline2D3D mRaster;    ///< a spline of the maximal size which represents the input function
  std::vector<float> mRasterData; ///< function values for the mRaster, with corrected edges

  SemiregularSpline2D3D mSpline;    ///< current spline
  std::vector<float> mSplineData; ///< function values for the mSpline, with corrected edges

  SemiregularSpline2D3D mTrySpline;    ///< spline to test the cost of the current action
  std::vector<float> mTrySplineData; ///< function values for the mTrySpline, with corrected edges
  
  int mCalibrationStage = 2; ///< current stage of the calibration: 0: V axis, 1: U axis, 2:stop
  float mMaxDeviation = 0.1; ///< maximal tolerated deviation between the spline and the input function
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
