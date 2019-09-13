// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  CompactSplineHelper.cxx
/// \brief Implementation of CompactSplineHelper class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE) && !defined(GPUCA_ALIROOT_LIB)

#include "CompactSplineHelper.h"
#include "TMath.h"
#include "TMatrixD.h"
#include <TVectorD.h>

using namespace GPUCA_NAMESPACE::gpu;

CompactSplineHelper::CompactSplineHelper() : mError() {}

std::unique_ptr<float[]> CompactSplineHelper::create(const CompactSplineIrregular1D& spline, std::function<float(float)> F, int nAxiliaryPoints)
{
  if (!spline.isConstructed()) {
    storeError(-1, "CompactSplineHelper::create: input spline is not constructed");
    return nullptr;
  }

  if (nAxiliaryPoints < 2) {
    nAxiliaryPoints = 2;
  }

  const int nKnots = spline.getNumberOfKnots();
  const int nPar = 2 * nKnots;
  std::unique_ptr<float[]> data(new float[nPar]);

  TMatrixDSym A(nPar);
  TVectorD b(nPar);

  A.Zero();
  b.Zero();

  auto addPoint = [&](double u, double f) {
    int i = spline.getKnotIndex(u);
    const CompactSplineIrregular1D::Knot& knot0 = spline.getKnot(i);
    const CompactSplineIrregular1D::Knot& knot1 = spline.getKnot(i + 1);
    double l = knot1.u - knot0.u;
    double x = (u - knot0.u) * knot0.Li; // scaled u
    double x2 = x * x;
    double xm1 = x - 1.;
    double cf1 = x2 * (3. - 2. * x);
    double cf0 = 1. - cf1;
    double cz0 = x * xm1 * xm1 * l;
    double cz1 = x2 * xm1 * l;

    i *= 2;
    A(i, i) += cf0 * cf0;
    A(i + 1, i) += cf0 * cz0;
    A(i + 2, i) += cf0 * cf1;
    A(i + 3, i) += cf0 * cz1;
    b[i] += cf0 * f;

    A(i + 1, i + 1) += cz0 * cz0;
    A(i + 2, i + 1) += cz0 * cf1;
    A(i + 3, i + 1) += cz0 * cz1;
    b[i + 1] += cz0 * f;

    A(i + 2, i + 2) += cf1 * cf1;
    A(i + 3, i + 2) += cf1 * cz1;
    b[i + 2] += cf1 * f;

    A(i + 3, i + 3) += cz1 * cz1;
    b[i + 3] += cz1 * f;
  };

  int nSteps = nAxiliaryPoints + 1;

  for (int i = 0; i < nKnots - 1; ++i) {
    const CompactSplineIrregular1D::Knot& knot0 = spline.getKnot(i);
    const CompactSplineIrregular1D::Knot& knot1 = spline.getKnot(i + 1);
    double u = knot0.u;
    double du = (knot1.u - u) / nSteps;
    for (int i = 0; i < nSteps; ++i, u += du) {
      addPoint(u, F(u));
    }
  }

  addPoint(1., F(1.)); // last knot

  // copy symmetric elements

  for (int i = 0; i < nPar; i++) {
    for (int j = i + 1; j < nPar; j++) {
      A(i, j) = A(j, i);
    }
  }

  TVectorD c = A.Invert() * b;
  for (int i = 0; i < nPar; i++) {
    data[i] = c[i];
  }

  return data;
}

std::unique_ptr<float[]> CompactSplineHelper::createClassical(const CompactSplineIrregular1D& spline, std::function<float(float)> F)
{
  if (!spline.isConstructed()) {
    storeError(-1, "CompactSplineHelper::create: input spline is not constructed");
    return nullptr;
  }

  const int nKnots = spline.getNumberOfKnots();
  std::unique_ptr<float[]> data(new float[2 * nKnots]);

  TMatrixD A(nKnots, nKnots);
  TVectorD b(nKnots);

  A.Zero();
  b.Zero();

  for (int i = 0; i < nKnots; ++i) {
    const CompactSplineIrregular1D::Knot& knot = spline.getKnot(i);
    double u = knot.u;
    data[2 * i] = F(u);
  }

  /*
    const CompactSplineIrregular1D::Knot& knot0 = spline.getKnot(i);
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
    const CompactSplineIrregular1D::Knot& knot0 = spline.getKnot(0);
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
    const CompactSplineIrregular1D::Knot& knot0 = spline.getKnot(nKnots - 2);
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
    const CompactSplineIrregular1D::Knot& knot0 = spline.getKnot(i - 1);
    double cf1 = (6 - 12) * knot0.Li * knot0.Li;
    double cz0 = (6 - 4) * knot0.Li;
    double cz1_0 = (6 - 2) * knot0.Li;
    // f''(u) = cf1*(f1-f0) + cz0*z0 + cz1*z1;

    const CompactSplineIrregular1D::Knot& knot1 = spline.getKnot(i);
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

#endif
