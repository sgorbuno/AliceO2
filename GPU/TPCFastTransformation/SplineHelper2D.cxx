// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  Spline182.cxx
/// \brief Implementation of SplineHelper2D class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE) && !defined(GPUCA_ALIROOT_LIB)

#include "SplineHelper2D.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompBK.h"

using namespace GPUCA_NAMESPACE::gpu;

SplineHelper2D::SplineHelper2D() : mError(), mSpline()
{
  mSpline.constructRegular(2, 2);
}

int SplineHelper2D::setSpline(const Spline2D& spline, int nAxiliaryPointsU, int nAxiliaryPointsV)
{
  // Prepare creation of 2D irregular spline
  // The should be at least one (better, two) axiliary measurements on each segnment between two knots and at least 2*nKnots measurements in total
  // Returns 0 when the spline can not be constructed with the given nAxiliaryPoints

  int ret = 0;

  if (!spline.isConstructed()) {
    ret = storeError(-1, "SplineHelper2D::setSpline2D: input spline is not constructed");
    mSpline.constructRegular(2, 2);
  } else {
    mSpline.cloneFromObject(spline, nullptr);
  }
  if (mHelperU.setSpline(mSpline.getGridU(), nAxiliaryPointsU) != 0) {
    ret = storeError(-2, "SplineHelper2D::setSpline2D: error by setting U axis");
  }

  if (mHelperV.setSpline(mSpline.getGridV(), nAxiliaryPointsV) != 0) {
    ret = storeError(-3, "SplineHelper2D::setSpline2D: error by setting V axis");
  }

  return ret;
}

void SplineHelper2D::constructParameters(int Ndim, const float F[/*getNumberOfMeasurements() x Ndim*/], float parameters[/*mSpline.getNumberOfParameters(Ndim) */]) const
{
  // Create 2D irregular spline in a compact way

  const int Ndim2 = 2 * Ndim;
  const int Ndim3 = 3 * Ndim;
  const int Ndim4 = 4 * Ndim;

  int nMeasurementsU = getNumberOfMeasurementsU();
  int nMeasurementsV = getNumberOfMeasurementsV();

  int nKnotsU = mSpline.getGridU().getNumberOfKnots();
  int nKnotsV = mSpline.getGridV().getNumberOfKnots();

  std::unique_ptr<float[]> rotF(new float[nMeasurementsU * nMeasurementsV * Ndim]); // U Measurements x V Measurements :  rotated F for one output dimension
  std::unique_ptr<float[]> Dv(new float[nKnotsV * nMeasurementsU * Ndim]);          // V knots x U Measurements

  std::unique_ptr<float[]> parU(new float[mHelperU.getSpline().getNumberOfParameters(Ndim)]);
  std::unique_ptr<float[]> parV(new float[mHelperV.getSpline().getNumberOfParameters(Ndim)]);

  // get the function values and U derivatives at knots from the U splines

  for (int ipu = 0; ipu < nMeasurementsU; ipu++) {
    for (int ipv = 0; ipv < nMeasurementsV; ipv++) {
      for (int dim = 0; dim < Ndim; dim++) {
        rotF[Ndim * (ipu * nMeasurementsV + ipv) + dim] = 0; //F[Ndim * (ipv * nMeasurementsU + ipu) + dim];
      }
    }
  }

  for (int iKnotV = 0; iKnotV < nKnotsV; ++iKnotV) {
    int ipv = mHelperV.getKnotMeasurement(iKnotV);
    const float* Frow = &(F[Ndim * ipv * nMeasurementsU]);
    mHelperU.constructParametersGradually(Ndim, Frow, parU.get());

    for (int iKnotU = 0; iKnotU < nKnotsU; iKnotU++) {
      for (int dim = 0; dim < Ndim; dim++) {
        parameters[Ndim4 * (iKnotV * nKnotsU + iKnotU) + dim] = parU[2 * iKnotU + dim];                // store f for all the knots
        parameters[Ndim4 * (iKnotV * nKnotsU + iKnotU) + Ndim2 + dim] = parU[2 * iKnotU + Ndim + dim]; // store f'u for all the knots
      }
    }

    // recalculate F values for all ipu Measurements at V = ipv
    for (int ipu = 0; ipu < nMeasurementsU; ipu++) {
      float splineF[Ndim];
      float u = mHelperU.getMeasurementPoint(ipu).u;
      mSpline.getGridU().interpolate(Ndim, parU.get(), u, splineF);
      for (int dim = 0; dim < Ndim; dim++) {
        rotF[(ipu * nMeasurementsV + ipv) * Ndim + dim] = splineF[dim];
      }
    }
  }

  for (int ipu = 0; ipu < nMeasurementsU; ipu++) {
    const float* Fcolumn = &(rotF[ipu * nMeasurementsV * Ndim]);
    mHelperV.constructParametersGradually(Ndim, Fcolumn, parV.get());
    for (int iKnotV = 0; iKnotV < nKnotsV; iKnotV++) {
      for (int dim = 0; dim < Ndim; dim++) {
        //int ipv = mHelperV.getKnotMeasurement(iKnotV);
        float dv = parV[(iKnotV * 2 + 1) * Ndim + dim];
        Dv[(iKnotV * nMeasurementsU + ipu) * Ndim + dim] = dv;
      }
    }
  }

  // fit F'v and f''_vu

  for (int iKnotV = 0; iKnotV < nKnotsV; ++iKnotV) {
    const float* measurements = &(Dv[iKnotV * nMeasurementsU * Ndim]);
    mHelperU.constructParametersGradually(Ndim, measurements, parU.get());
    for (int iKnotU = 0; iKnotU < nKnotsU; iKnotU++) {
      for (int dim = 0; dim < Ndim; dim++) {
        parameters[iKnotV * Ndim4 * nKnotsU + iKnotU * Ndim4 + Ndim + dim] = parU[2 * iKnotU + 0];  // store f'v for all the knots
        parameters[iKnotV * Ndim4 * nKnotsU + iKnotU * Ndim4 + Ndim3 + dim] = parU[2 * iKnotU + 1]; // store f''vu for all the knots
      }
    }
  }
}

  /*

    std::unique_ptr<float[]> SplineHelper2D::create(const Spline2D& spline, std::function<void(float, float, float&, float&, float&)> F, int nAxiliaryPoints)
    {
      if (!spline.isConstructed()) {
        storeError(-1, "SplineHelper2D::create: input spline is not constructed");
        return nullptr;
      }

      if (nAxiliaryPoints < 2) {
        nAxiliaryPoints = 2;
      }

      const int nKnots = spline.getNumberOfKnots();
      const int nPar = 4 * nKnots;
      std::unique_ptr<float[]> parameters(new float[nPar]);
      for (int i = 0; i < nPar; i++)
        parameters[i] = 0.f;
      return parameters;
    }
*/

#endif
