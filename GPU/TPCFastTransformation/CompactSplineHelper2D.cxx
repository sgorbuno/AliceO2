// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  CompactSpline182.cxx
/// \brief Implementation of CompactSplineHelper2D class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE) && !defined(GPUCA_ALIROOT_LIB)

#include "CompactSplineHelper2D.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompBK.h"

using namespace GPUCA_NAMESPACE::gpu;

CompactSplineHelper2D::CompactSplineHelper2D() : mError() {}

int CompactSplineHelper2D::setSpline(const CompactSpline2D& spline, int nAxiliaryPointsU, int nAxiliaryPointsV)
{
  // Prepare creation of 2D irregular spline
  // The should be at least one (better, two) axiliary data point on each segnment between two knots and at least 2*nKnots data points in total
  // Returns 0 when the spline can not be constructed with the given nAxiliaryPoints

  int ret = 0;

  if (!spline.isConstructed()) {
    ret = storeError(-1, "CompactSplineHelper2D::setSpline2D: input spline is not constructed");
  }

  mSpline.cloneFromObject(spline, nullptr);

  if (mHelperU.setSpline(spline.getGridU(), nAxiliaryPointsU) != 0) {
    ret = storeError(-2, "CompactSplineHelper2D::setSpline2D: error by setting U axis");
  }

  if (mHelperV.setSpline(spline.getGridV(), nAxiliaryPointsV) != 0) {
    ret = storeError(-3, "CompactSplineHelper2D::setSpline2D: error by setting V axis");
  }

  return ret;
}

void CompactSplineHelper2D::constructSpline(const float inF[/*getNdataPoints()*/], float outSplineData[/*getNparameters()*/]) const
{
  // Create 2D irregular spline in a compact way

  int nPointsU = getNdataPointsU();
  int nPointsV = getNdataPointsV();

  int nKnotsU = mSpline.getGridU().getNumberOfKnots();
  int nKnotsV = mSpline.getGridV().getNumberOfKnots();

  float mapF[nPointsU * nPointsV]; // V points x U points :  rotated inF for one dimension
  float mapFv[nKnotsV * nPointsU]; // U points x V knots

  float pointsU[nPointsU];

  float dataU[mHelperU.getNparameters()];
  float dataV[mHelperV.getNparameters()];

  for (int dim = 0; dim < 3; dim++) { // loop over dimensions

    // get the function values and U derivatives at knots from the U splines

    for (int ipu = 0; ipu < nPointsU; ipu++) {
      for (int ipv = 0; ipv < nPointsV; ipv++) {
        mapF[ipu * nPointsV + ipv] = inF[3 * (ipv * nPointsU + ipu) + dim];
      }
    }

    for (int iKnotV = 0; iKnotV < nKnotsV; ++iKnotV) {
      int ipv = mHelperV.getKnotPoint(iKnotV);
      const float* inFrow = &(inF[ipv * 3 * nPointsU]);
      for (int ipu = 0; ipu < nPointsU; ipu++) {
        pointsU[ipu] = inFrow[3 * ipu + dim];
      }
      mHelperU.constructDataGradually(1, pointsU, dataU);

      for (int iKnotU = 0; iKnotU < nKnotsU; iKnotU++) {
        outSplineData[iKnotV * 12 * nKnotsU + iKnotU * 12 + dim] = dataU[2 * iKnotU + 0];     // store f for all the knots
        outSplineData[iKnotV * 12 * nKnotsU + iKnotU * 12 + 6 + dim] = dataU[2 * iKnotU + 1]; // store f'u for all the knots
      }

      // recalculate F values for all ipu points at V = ipv
      for (int ipu = 0; ipu < nPointsU; ipu++) {
        float splineF;
        float u = ipu * nKnotsU / nPointsU;
        mSpline.getGridU().getSpline<1>(dataU, u, &splineF);
        mapF[ipu * nPointsV + ipv] = splineF;
      }
    }

    for (int ipu = 0; ipu < nPointsU; ipu++) {
      float* points = &(mapF[ipu * nPointsV]);
      mHelperV.constructDataGradually(1, points, dataV);
      for (int iKnotV = 0; iKnotV < nKnotsV; iKnotV++) {
        int ipv = mHelperV.getKnotPoint(iKnotV);
        float dv = dataV[iKnotV * 2 + 1];
        mapFv[iKnotV * nPointsU + ipu] = dv;
      }
    }

    // fit F'v and f''_vu

    for (int iKnotV = 0; iKnotV < nKnotsV; ++iKnotV) {
      float* points = &(mapFv[iKnotV * nPointsU]);
      mHelperU.constructDataGradually(1, pointsU, dataU);
      for (int iKnotU = 0; iKnotU < nKnotsU; iKnotU++) {
        outSplineData[iKnotV * 12 * nKnotsU + iKnotU * 12 + 3 + dim] = dataU[2 * iKnotU + 0]; // store f'v for all the knots
        outSplineData[iKnotV * 12 * nKnotsU + iKnotU * 12 + 9 + dim] = dataU[2 * iKnotU + 1]; // store f''vu for all the knots
      }
    }

  } // dimensions
}

/*

    std::unique_ptr<float[]> CompactSplineHelper2D::create(const CompactSpline2D& spline, std::function<void(float, float, float&, float&, float&)> F, int nAxiliaryPoints)
    {
      if (!spline.isConstructed()) {
        storeError(-1, "CompactSplineHelper2D::create: input spline is not constructed");
        return nullptr;
      }

      if (nAxiliaryPoints < 2) {
        nAxiliaryPoints = 2;
      }

      const int nKnots = spline.getNumberOfKnots();
      const int nPar = 4 * nKnots;
      std::unique_ptr<float[]> data(new float[nPar]);
      for (int i = 0; i < nPar; i++)
        data[i] = 0.f;
      return data;
    }
*/

#endif
