// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file SemiregularSpline2D3DCalibrator.cxx
/// \brief Implementation of SemiregularSpline2D3DCalibrator class
///
/// \author  Anna Bartsch
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>


#include "SemiregularSpline2D3D.h"
#include "SemiregularSpline2D3DCalibrator.h"


namespace GPUCA_NAMESPACE
{
namespace gpu
{

SemiregularSpline2D3DCalibrator::SemiregularSpline2D3DCalibrator()
{
  /// Default constructor
  setRasterSize(5, 5);
  setMaxNKnots(5, 5);
}

void SemiregularSpline2D3DCalibrator::setRasterSize(int nTicksU, int nTicksV)
{
  /// set maximal size of the spline grid

  if (nTicksU < MinNKnots)
    nTicksU = MinNKnots;  

  if (nTicksV <MinNKnots) 
    nTicksV = MinNKnots;
  
  int nTicksPerRow[nTicksV];
  for( int i=0; i<nTicksV; i++ ) nTicksPerRow[i] = nTicksU;
  mRaster.construct(nTicksV, nTicksPerRow );
}

void SemiregularSpline2D3DCalibrator::setMaxNKnots(int nKnotsU, int nKnotsV)
{
  /// set maximal size of the spline grid

  mMaxNKnots[0] = nKnotsU; //umbenennen?
  mMaxNKnots[1] = nKnotsV;

  for (int uv = 0; uv < 2; uv++) {
    if (mMaxNKnots[uv] < MinNKnots)
      mMaxNKnots[uv] = MinNKnots;
  }
}

void SemiregularSpline2D3DCalibrator::startCalibration(void (*F)(float, float, float&, float&, float&))
{
  // initialize everything for the calibration

  // fill the raster data
  mRasterData.resize(mRaster.getNumberOfKnots() * 3);

  for (int i = 0; i < mRaster.getNumberOfKnots(); i++) {
    float u = 0, v = 0, fx = 0, fy = 0, fz = 0;
    mRaster.getKnotUV(i, u, v);
    F(u, v, fx, fy, fz);
    mRasterData[3 * i + 0] = fx;
    mRasterData[3 * i + 1] = fy;
    mRasterData[3 * i + 2] = fz;
  }

  mRaster.correctEdges(mRasterData.data());

  mNKnotsV = mMaxNKnots[1];
  mNKnotsU.resize( mNKnotsV );

  for (int i=0; i<mNKnotsV; i++) mNKnotsU[i] = mMaxNKnots[0];

  createCurrentSpline();
  mCalibrationStage = 0;
}

void SemiregularSpline2D3DCalibrator::createCurrentSpline()
{
  createSpline(mSpline, mSplineData);
}

void SemiregularSpline2D3DCalibrator::createTrySpline()
{
  createSpline(mTrySpline, mTrySplineData);
}

void SemiregularSpline2D3DCalibrator::createSpline(SemiregularSpline2D3D& sp, std::vector<float>& data)
{
  // recreate a spline with  respect to knots

  sp.construct(mNKnotsV, mNKnotsU.data() );  

  data.resize(sp.getNumberOfKnots() * 3);
  for (int i = 0; i < sp.getNumberOfKnots(); i++) {
    float u = 0, v = 0, fx = 0, fy = 0, fz = 0;
    sp.getKnotUV(i, u, v);
    mRaster.getSplineVec(mRasterData.data(), u, v, fx, fy, fz);
    data[3 * i + 0] = fx;
    data[3 * i + 1] = fy;
    data[3 * i + 2] = fz;
  }
  sp.correctEdges(data.data());
}


 bool SemiregularSpline2D3DCalibrator::doCalibrationStep()
{
  // perform one step of the calibration

std::cout<<"calibration stage "<< mCalibrationStage<<" mNKnotsV= "<<mNKnotsV<<std::endl;
  if( mCalibrationStage == 0 ){ // try to adjust V axis    
    if( mNKnotsV > 5 ){
      mNKnotsV--;
      createTrySpline();
      double tryCost=0, tryMaxDev = 0;
      getCost(mTrySpline, mTrySplineData, tryCost, tryMaxDev );
      tryMaxDev = sqrt(tryMaxDev/3.);
      std::cout<<" try to decrease V: max "<<tryMaxDev<<", cost "<<tryCost<<std::endl; 
      if( tryMaxDev <= mMaxDeviation ) return 1;
      mNKnotsV++;
    }
    mCalibrationStage = 1;    
  }
  
  if( mCalibrationStage == 1 ){ // get the row with the best cost and remove a knot
    //createCurrentSpline(); 
    //double currentCost=0, currentMaxDev = 0;
    //getCost(mSpline, mSplineData, currentCost, currentMaxDev );    
    //if( currentMaxDev > mMaxDeviation ){
      //mCalibrationStage = 2;    
      //return 0;
    //}
    double bestCost=1.e100;
    double bestDev=0.;
    int bestRow=-1;
    for( int i=0; i<mNKnotsV; i++ ){
      if( mNKnotsU[i]>= MinNKnots ) continue;
      mNKnotsU[i]--;
      createTrySpline();
      mNKnotsU[i]++;
      double tryCost=0, tryMaxDev=0;
      getCost(mTrySpline, mTrySplineData, tryCost, tryMaxDev );
      tryMaxDev = sqrt(tryMaxDev/3.);
	std::cout<<" U row "<<i<<": max "<<tryMaxDev<<", cost "<<tryCost<<std::endl; 
      if (tryMaxDev < mMaxDeviation && tryCost < bestCost){
        bestCost = tryCost;
	bestDev= tryMaxDev;
        bestRow = i;
      }
    }
    if (bestRow >= 0) {
      mNKnotsU[bestRow]--;
      createCurrentSpline();
      std::cout << "bestRow: " << bestRow << " bestCost: " << bestCost << " bestDev: " << bestDev << std::endl;
      return 1;
    }
  }
  return 0;
}


void SemiregularSpline2D3DCalibrator::getCost(const SemiregularSpline2D3D& spline, const std::vector<float>& data, double &cost, double &maxDeviation) const
{
  // get the maximal Deviation and the sum cost

  maxDeviation=0.;
  cost=0.;
  for ( int i=0; i < mRaster.getNumberOfKnots(); i++){
    float u=0., v=0.;
    mRaster.getKnotUV(i, u, v);     
    float fx0, fy0, fz0, fx, fy, fz;
    mRaster.getSplineVec(mRasterData.data(), u, v, fx0,fy0,fy0);
    spline.getSplineVec(data.data(), u, v, fx, fy, fz);
    double dx = fx-fx0;
    double dy = fy-fy0;
    double dz = fz-fz0;
    double d = dx*dx + dy*dy + dz*dz;
    if (d > maxDeviation) maxDeviation = d;
    cost+=d;
  }
}

std::unique_ptr<float[]> SemiregularSpline2D3DCalibrator::calibrateSpline(SemiregularSpline2D3D& spline_uv, void (*F)(float, float, float&, float&, float&))
{
  // main method: spline calibration

  startCalibration(F);
  while (doCalibrationStep());
  createCurrentSpline();
  spline_uv.cloneFromObject(mSpline, nullptr);
  std::unique_ptr<float[]> tmp(new float[mSpline.getNumberOfKnots()]);
  for (int i = 0; i < mSpline.getNumberOfKnots(); i++) {
    tmp[i] = mSplineData[i];
  }
  return tmp;
  }

} // namespace gpu
} // namespace GPUCA_NAMESPACE
