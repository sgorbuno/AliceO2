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
/// \brief Implementation of CompactSplineHelper1D class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE) && !defined(GPUCA_ALIROOT_LIB)

#include "CompactSplineHelper1D.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompBK.h"

using namespace GPUCA_NAMESPACE::gpu;

CompactSplineHelper1D::CompactSplineHelper1D() : mError() {}

std::unique_ptr<float[]> CompactSplineHelper1D::constructDataClassical1D(const CompactSpline1D& spline, std::function<float(float)> F, float uMin, float uMax)
{
  // Create 1D->1D spline in a classical way:
  // set slopes at the knots such, that the second derivative of the spline stays continious.
  //
  if (!spline.isConstructed()) {
    storeError(-1, "CompactSplineHelper1D::create: input spline is not constructed");
    return nullptr;
  }

  const int nKnots = spline.getNumberOfKnots();
  std::unique_ptr<float[]> data(new float[2 * nKnots]);

  TMatrixD A(nKnots, nKnots);
  TVectorD b(nKnots);

  A.Zero();
  b.Zero();
  double scale = (uMax - uMin) / ((double)nKnots - 1.);
  for (int i = 0; i < nKnots; ++i) {
    const CompactSpline1D::Knot& knot = spline.getKnot(i);
    double u = knot.u;
    data[2 * i] = F(uMin + u * scale);
  }

  /*
    const CompactSpline1D::Knot& knot0 = spline.getKnot(i);
    double x = (u - knot0.u) * knot0.Li; // scaled u    
    double cf1 = (6 - 12*x)*knot0.Li*knot0.Li;
    double cz0 = (6*x-4)*knot0.Li;
    double cz1 = (6*x-2)*knot0.Li;
    // f''(u) = cf1*(f1-f0) + cz0*z0 + cz1*z1;
   */

  // second derivative at knot0 is 0
  {
    double f0 = data[2 * 0];
    double f1 = data[2 * 1];
    const CompactSpline1D::Knot& knot0 = spline.getKnot(0);
    double cf1 = (6) * knot0.Li * knot0.Li;
    double cz0 = (-4) * knot0.Li;
    double cz1 = (-2) * knot0.Li;
    // f''(u) = cf1*(f1-f0) + cz0*z0 + cz1*z1;
    A(0, 0) = cz0;
    A(0, 1) = cz1;
    b(0) = -cf1 * (f1 - f0);
  }

  // second derivative at knot nKnots-1  is 0
  {
    double f0 = data[2 * (nKnots - 2)];
    double f1 = data[2 * (nKnots - 1)];
    const CompactSpline1D::Knot& knot0 = spline.getKnot(nKnots - 2);
    double cf1 = (6 - 12) * knot0.Li * knot0.Li;
    double cz0 = (6 - 4) * knot0.Li;
    double cz1 = (6 - 2) * knot0.Li;
    // f''(u) = cf1*(f1-f0) + cz0*z0 + cz1*z1;
    A(nKnots - 1, nKnots - 2) = cz0;
    A(nKnots - 1, nKnots - 1) = cz1;
    b(nKnots - 1) = -cf1 * (f1 - f0);
  }

  // second derivative at other knots is same from the left and from the right
  for (int i = 1; i < nKnots - 1; i++) {
    double f0 = data[2 * (i - 1)];
    double f1 = data[2 * (i)];
    double f2 = data[2 * (i + 1)];
    const CompactSpline1D::Knot& knot0 = spline.getKnot(i - 1);
    double cf1 = (6 - 12) * knot0.Li * knot0.Li;
    double cz0 = (6 - 4) * knot0.Li;
    double cz1_0 = (6 - 2) * knot0.Li;
    // f''(u) = cf1*(f1-f0) + cz0*z0 + cz1*z1;

    const CompactSpline1D::Knot& knot1 = spline.getKnot(i);
    double cf2 = (6) * knot1.Li * knot1.Li;
    double cz1_1 = (-4) * knot1.Li;
    double cz2 = (-2) * knot1.Li;
    // f''(u) = cf2*(f2-f1) + cz1_1*z1 + cz2*z2;
    A(i, i - 1) = cz0;
    A(i, i) = cz1_0 - cz1_1;
    A(i, i + 1) = -cz2;
    b(i) = -cf1 * (f1 - f0) + cf2 * (f2 - f1);
  }

  TVectorD c = A.Invert() * b;
  for (int i = 0; i < nKnots; i++) {
    data[2 * i + 1] = c[i];
  }

  return data;
}

int CompactSplineHelper1D::setSpline(const CompactSpline1D& spline, int nAxiliaryPoints)
{
  // Prepare creation of 1D irregular spline in a compact way:
  // fit all the spline parameters (which are the spline values and the slopes at the knots) to multiple data points.
  // The should be at least one (better, two) axiliary data point on each segnment between two knots and at least 2*nKnots data points in total
  // Returns 0 when the spline can not be constructed with the given nAxiliaryPoints

  int ret = 0;

  mNKnots = spline.getNumberOfKnots();
  int nPoints = 0;
  if (!spline.isConstructed()) {
    ret = storeError(-1, "CompactSplineHelper1D::setSpline: input spline is not constructed");
    mNKnots = 2;
    nAxiliaryPoints = 2;
    nPoints = 4;
  } else {
    if (nAxiliaryPoints < 1) {
      ret = storeError(-2, "CompactSplineHelper1D::setSpline: nAxiliaryPoints<1, increase to 1 ");
      nAxiliaryPoints = 1;
    }

    nPoints = 1 + spline.getUmax() + spline.getUmax() * nAxiliaryPoints;

    if (nPoints < 2 * spline.getNumberOfKnots()) {
      nAxiliaryPoints = 2;
      nPoints = 1 + spline.getUmax() + spline.getUmax() * nAxiliaryPoints;
      ret = storeError(-3, "CompactSplineHelper1D::setSpline: too few nAxiliaryPoints, increase to 2");
    }
  }

  const int nPar = getNparameters();

  mPoints.resize(nPoints);
  mKnotPoints.resize(mNKnots);
  for (int i = 0; i < mNKnots; ++i) {
    const CompactSpline1D::Knot& knot = spline.getKnot(i);
    int iu = (int)(knot.u + 0.1f);
    mKnotPoints[i] = iu * (1 + nAxiliaryPoints);
  }

  TMatrixDSym A(nPar);
  A.Zero();

  double scalePoints2Knots = ((double)spline.getUmax()) / (nPoints - 1.);
  for (int i = 0; i < nPoints; ++i) {
    Point& p = mPoints[i];
    double u = i * scalePoints2Knots;
    int iKnot = spline.getKnotIndex(u);
    const CompactSpline1D::Knot& knot0 = spline.getKnot(iKnot);
    const CompactSpline1D::Knot& knot1 = spline.getKnot(iKnot + 1);
    double l = knot1.u - knot0.u;
    double x = (u - knot0.u) * knot0.Li; // scaled u
    double x2 = x * x;
    double xm1 = x - 1.;

    p.iKnot = iKnot;
    p.cf1 = x2 * (3. - 2. * x);
    p.cf0 = 1. - p.cf1;
    p.cz0 = x * xm1 * xm1 * l;
    p.cz1 = x2 * xm1 * l;

    int j = iKnot * 2;
    A(j + 0, j + 0) += p.cf0 * p.cf0;
    A(j + 1, j + 0) += p.cf0 * p.cz0;
    A(j + 2, j + 0) += p.cf0 * p.cf1;
    A(j + 3, j + 0) += p.cf0 * p.cz1;

    A(j + 1, j + 1) += p.cz0 * p.cz0;
    A(j + 2, j + 1) += p.cz0 * p.cf1;
    A(j + 3, j + 1) += p.cz0 * p.cz1;

    A(j + 2, j + 2) += p.cf1 * p.cf1;
    A(j + 3, j + 2) += p.cf1 * p.cz1;

    A(j + 3, j + 3) += p.cz1 * p.cz1;
  }

  // copy symmetric matrix elements

  for (int i = 0; i < nPar; i++) {
    for (int j = i + 1; j < nPar; j++) {
      A(i, j) = A(j, i);
    }
  }

  TMatrixDSym Z(mNKnots);
  mMatrixFastF.resize(mNKnots * mNKnots);
  for (int i = 0, k = 0; i < mNKnots; i++) {
    for (int j = 0; j < mNKnots; j++, k++) {
      mMatrixFastF[k] = A(i * 2 + 1, j * 2);
      Z(i, j) = A(i * 2 + 1, j * 2 + 1);
    }
  }

  {
    TDecompBK bk(A, 0);
    bool ok = bk.Invert(A);

    if (!ok) {
      ret = storeError(-4, "CompactSplineHelper1D::setSpline: internal error - can not invert the matrix");
      A.Zero();
    }
    mMatrixI.resize(nPar * nPar);
    for (int i = 0, k = 0; i < nPar; i++) {
      for (int j = 0; j < nPar; j++, k++) {
        mMatrixI[k] = A(i, j);
      }
    }
  }

  {
    TDecompBK bk(Z, 0);
    if (!bk.Invert(Z)) {
      ret = storeError(-5, "CompactSplineHelper1D::setSpline: internal error - can not invert the matrix");
      Z.Zero();
    }
    mMatrixFastI.resize(mNKnots * mNKnots);
    for (int i = 0, k = 0; i < mNKnots; i++) {
      for (int j = 0; j < mNKnots; j++, k++) {
        mMatrixFastI[k] = Z(i, j);
      }
    }
  }

  return ret;
}

void CompactSplineHelper1D::constructData1D( const float inF[/*N Data Points*/], float outSplineData[/*N Spline Parameters*/]) const
{
  // Create 1D irregular spline in a compact way

  const int nPar = getNparameters();

  double b[nPar];
  for (int i = 0; i < nPar; i++)
    b[i] = 0.;

  for (int i = 0; i < getNdataPoints(); ++i) {
    const Point& p = mPoints[i];
    double* bb = &(b[p.iKnot * 2]);
    double f = (double)inF[i];
    bb[0] += f * p.cf0;
    bb[1] += f * p.cz0;
    bb[2] += f * p.cf1;
    bb[3] += f * p.cz1;
  }

  const double* row = mMatrixI.data();

  for (int i = 0; i < nPar; i++, row += nPar) {
    double s = 0.;
    for (int j = 0; j < nPar; j++) {
      s += row[j] * b[j];
    }
    outSplineData[i] = (float)s;
  }
}

std::unique_ptr<float[]> CompactSplineHelper1D::constructData1D(const CompactSpline1D& spline, std::function<float(float)> F, float uMin, float uMax, int nAxiliaryPoints)
{
  // Create 1D spline in a compact way for the input function F.
  // nAxiliaryPoints: number of data points between the spline knots (should be at least 2)

  if (!spline.isConstructed()) {
    storeError(-1, "CompactSplineHelper1D::constructData: input spline is not constructed");
    return nullptr;
  }
  if (nAxiliaryPoints < 2) {
    nAxiliaryPoints = 2;
  }
  int err = setSpline(spline, nAxiliaryPoints);
  if (err != 0)
    return nullptr;

  std::vector<float> vF(getNdataPoints());

  double scale = (uMax - uMin) / ((double)getNdataPoints() - 1.);
  for (int i = 0; i < getNdataPoints(); i++) {
    vF[i] = F(uMin + i * scale);
  }

  std::unique_ptr<float[]> splineData(new float[getNparameters()]);
  constructData1D(vF.data(), splineData.get());
  return splineData;
}

void CompactSplineHelper1D::constructDataGradually(int Ndim, const float inF[/*N Data Points x Ndim */], float outSplineData[/*N Spline Parameters * Ndim*/]) const
{
  // Create 1D irregular spline in a compact way

  for (int i = 0; i < mNKnots; ++i) { // set F values at knots
    int ip = mKnotPoints[i];
    for (int d = 0; d < Ndim; d++) {
      outSplineData[2 * i * Ndim + d] = inF[ip * Ndim + d];
    }
  }

  double b[mNKnots * Ndim];
  for (int i = 0; i < mNKnots * Ndim; i++) {
    b[i] = 0.;
  }

  for (int i = 0; i < getNdataPoints(); ++i) {
    const Point& p = mPoints[i];
    for (int d = 0; d < Ndim; d++) {
      double f = (double)inF[i * Ndim + d];
      b[(p.iKnot + 0) * Ndim + d] += f * p.cz0;
      b[(p.iKnot + 1) * Ndim + d] += f * p.cz1;
    }
  }

  const double* row = mMatrixFastF.data();
  for (int i = 0; i < mNKnots; ++i, row += mNKnots) {
    double s[Ndim];
    for (int d = 0; d < Ndim; d++)
      s[d] = 0.;
    for (int j = 0; j < mNKnots; ++j) {
      for (int d = 0; d < Ndim; d++)
        s[d] += row[j] * outSplineData[2 * j * Ndim + d];
    }
    for (int d = 0; d < Ndim; d++)
      b[i * Ndim + d] -= s[d];
  }

  row = mMatrixFastI.data();
  for (int i = 0; i < mNKnots; ++i, row += mNKnots) {
    double s[Ndim];
    for (int d = 0; d < Ndim; d++) {
      s[d] = 0.;
    }
    for (int j = 0; j < mNKnots; ++j) {
      for (int d = 0; d < Ndim; d++) {
        s[d] += row[j] * b[j * Ndim + d];
      }
    }
    for (int d = 0; d < Ndim; d++) {
      outSplineData[2 * (i + 1) * Ndim + d] = (float)s[d];
    }
  }
}


#endif
