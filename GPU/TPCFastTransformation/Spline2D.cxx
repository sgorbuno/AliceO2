// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  Spline2D.cxx
/// \brief Implementation of Spline2D class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#include "Spline2D.h"

#if !defined(GPUCA_GPUCODE)
#include <iostream>
#endif

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE) // code invisible on GPU and in the standalone compilation
#include "TRandom.h"
#include "Riostream.h"
#include "TMath.h"
#include "SplineHelper2D.h"
#include "TCanvas.h"
#include "TNtuple.h"
#endif

using namespace GPUCA_NAMESPACE::gpu;

Spline2D::Spline2D() : FlatObject(), mGridU(2), mGridV(2)
{
  /// Default constructor. Creates an empty uninitialised object
}

void Spline2D::destroy()
{
  /// See FlatObject for description
  mGridU.destroy();
  mGridV.destroy();
  FlatObject::destroy();
}

void Spline2D::cloneFromObject(const Spline2D& obj, char* newFlatBufferPtr)
{
  /// See FlatObject for description

  const char* oldFlatBufferPtr = obj.mFlatBufferPtr;

  FlatObject::cloneFromObject(obj, newFlatBufferPtr);

  char* bufferU = FlatObject::relocatePointer(oldFlatBufferPtr, mFlatBufferPtr, obj.mGridU.getFlatBufferPtr());
  char* bufferV = FlatObject::relocatePointer(oldFlatBufferPtr, mFlatBufferPtr, obj.mGridV.getFlatBufferPtr());

  mGridU.cloneFromObject(obj.mGridU, bufferU);
  mGridV.cloneFromObject(obj.mGridV, bufferV);
}

void Spline2D::moveBufferTo(char* newFlatBufferPtr)
{
/// See FlatObject for description
#ifndef GPUCA_GPUCODE
  char* oldFlatBufferPtr = mFlatBufferPtr;
  FlatObject::moveBufferTo(newFlatBufferPtr);
  char* currFlatBufferPtr = mFlatBufferPtr;
  mFlatBufferPtr = oldFlatBufferPtr;
  setActualBufferAddress(currFlatBufferPtr);
#endif
}

void Spline2D::setActualBufferAddress(char* actualFlatBufferPtr)
{
  /// See FlatObject for description
  char* bufferU = FlatObject::relocatePointer(mFlatBufferPtr, actualFlatBufferPtr, mGridU.getFlatBufferPtr());
  char* bufferV = FlatObject::relocatePointer(mFlatBufferPtr, actualFlatBufferPtr, mGridV.getFlatBufferPtr());
  mGridU.setActualBufferAddress(bufferU);
  mGridV.setActualBufferAddress(bufferV);
  FlatObject::setActualBufferAddress(actualFlatBufferPtr);
}

void Spline2D::setFutureBufferAddress(char* futureFlatBufferPtr)
{
  /// See FlatObject for description
  char* bufferU = relocatePointer(mFlatBufferPtr, futureFlatBufferPtr, mGridU.getFlatBufferPtr());
  char* bufferV = relocatePointer(mFlatBufferPtr, futureFlatBufferPtr, mGridV.getFlatBufferPtr());
  mGridU.setFutureBufferAddress(bufferU);
  mGridV.setFutureBufferAddress(bufferV);
  FlatObject::setFutureBufferAddress(futureFlatBufferPtr);
}

void Spline2D::construct(int numberOfKnotsU, const int knotsU[], int numberOfKnotsV, const int knotsV[])
{
  /// Constructor
  ///
  /// Number of created knots may differ from the input values:
  /// - Edge knots {0} and {Umax/Vmax} will be added if they are not present.
  /// - Duplicated knots, knots with a negative coordinate will be deleted
  /// - At least 2 knots for each axis will be created
  ///
  /// \param numberOfKnotsU     Number of knots in knotsU[] array
  /// \param knotsU             Array of knot positions (integer values)
  ///
  /// \param numberOfKnotsV     Number of knots in knotsV[] array
  /// \param knotsV             Array of knot positions (integer values)
  ///

  FlatObject::startConstruction();

  mGridU.construct(numberOfKnotsU, knotsU);
  mGridV.construct(numberOfKnotsV, knotsV);

  size_t vOffset = alignSize(mGridU.getFlatBufferSize(), mGridV.getBufferAlignmentBytes());

  FlatObject::finishConstruction(vOffset + mGridV.getFlatBufferSize());

  mGridU.moveBufferTo(mFlatBufferPtr);
  mGridV.moveBufferTo(mFlatBufferPtr + vOffset);
}

void Spline2D::constructRegular(int numberOfKnotsU, int numberOfKnotsV)
{
  /// Constructor for a regular spline
  /// \param numberOfKnotsU     U axis: Number of knots in knots[] array
  /// \param numberOfKnotsV     V axis: Number of knots in knots[] array
  ///

  FlatObject::startConstruction();

  mGridU.constructRegular(numberOfKnotsU);
  mGridV.constructRegular(numberOfKnotsV);

  size_t vOffset = alignSize(mGridU.getFlatBufferSize(), mGridV.getBufferAlignmentBytes());

  FlatObject::finishConstruction(vOffset + mGridV.getFlatBufferSize());

  mGridU.moveBufferTo(mFlatBufferPtr);
  mGridV.moveBufferTo(mFlatBufferPtr + vOffset);
}

void Spline2D::print() const
{
#if !defined(GPUCA_GPUCODE)
  std::cout << " Irregular Spline 2D: " << std::endl;
  std::cout << " grid U: " << std::endl;
  mGridU.print();
  std::cout << " grid V: " << std::endl;
  mGridV.print();
#endif
}

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE) // code invisible on GPU and in the standalone compilation

int Spline2D::test(bool draw)
{
  using namespace std;

  const int funcN = 2;
  double funcC[4 * funcN * funcN + 4];

  int nKnots = 4;
  const int nAxiliaryPoints = 5;
  int uMax = nKnots * 3;

  auto F = [&](float u, float v) -> float {
    double uu = u * TMath::Pi() / uMax;
    double vv = 0;//v * TMath::Pi() / uMax;
    double f = 0; //funcC[0]/2;
    for (int i = 1; i <= funcN; i++) {
      for (int j = 1; j <= funcN; j++) {
        f += funcC[4 * (i * funcN + j) + 0] * TMath::Cos(i * uu) * TMath::Cos(j * vv);
        f += funcC[4 * (i * funcN + j) + 1] * TMath::Cos(i * uu) * TMath::Sin(j * vv);
        f += funcC[4 * (i * funcN + j) + 2] * TMath::Sin(i * uu) * TMath::Cos(j * vv);
        f += funcC[4 * (i * funcN + j) + 3] * TMath::Sin(i * uu) * TMath::Sin(j * vv);
      }
    }
    return f;
  };

  TCanvas* canv = nullptr;
  TNtuple* nt = nullptr;
  TNtuple* knots = nullptr;

  auto ask = [&]() -> bool {
    if (!canv)
      return 0;
    canv->Update();
    cout << "type 'q ' to exit" << endl;
    std::string str;
    std::getline(std::cin, str);
    return (str != "q" && str != ".q");
  };

  std::cout
    << "Test 2D interpolation with the compact spline" << std::endl;

  int nTries = 100;

  if (draw) {
    canv = new TCanvas("cQA", "Spline1D  QA", 2000, 1000);
    nTries = 10000;
  }

  double statDf = 0;
  double statN = 0;

  for (int seed = 1; seed < nTries; seed++) {
    cout << "next try.." << endl;

    gRandom->SetSeed(seed);

    for (int i = 0; i < 4 * funcN * funcN + 4; i++) {
      funcC[i] = gRandom->Uniform(-1, 1);
    }

    SplineHelper2D helper;
    Spline2D spline;

    do {
      int knotsU[nKnots], knotsV[nKnots];
      knotsU[0] = 0;
      knotsV[0] = 0;
      double du = 1. * uMax / (nKnots - 1);
      for (int i = 1; i < nKnots; i++) {
        knotsU[i] = (int)(i * du); // + gRandom->Uniform(-du / 3, du / 3);
        knotsV[i] = (int)(i * du); // + gRandom->Uniform(-du / 3, du / 3);
      }
      knotsU[nKnots - 1] = uMax;
      knotsV[nKnots - 1] = uMax;
      spline.construct(nKnots, knotsU, nKnots, knotsV);

      if (nKnots != spline.getGridU().getNumberOfKnots() || nKnots != spline.getGridV().getNumberOfKnots()) {
        cout << "warning: n knots changed during the initialisation " << nKnots << " -> " << spline.getNumberOfKnots() << std::endl;
        continue;
      }
    } while (0);

    std::string err = FlatObject::stressTest(spline);
    if (!err.empty()) {
      cout << "error at FlatObject functionality: " << err << endl;
      return -1;
    } else {
      //cout << "flat object functionality is ok" << endl;
    }

    helper.setSpline(spline, nAxiliaryPoints, nAxiliaryPoints);

    std::unique_ptr<float[]> parameters(new float[spline.getNumberOfParameters(1)]);
    std::unique_ptr<float[]> mapF(new float[helper.getNumberOfMeasurements()]);

    int nPointsU = helper.getNumberOfMeasurementsU();
    int nPointsV = helper.getNumberOfMeasurementsV();
    for (int ipu = 0; ipu < nPointsU; ipu++) {
      for (int ipv = 0; ipv < nPointsV; ipv++) {
        mapF[ipu * nPointsV + ipv] = F(helper.getHelperU().getMeasurementPoint(ipu).u, helper.getHelperV().getMeasurementPoint(ipv).u);
      }
    }

    helper.constructParameters(1, mapF.get(), parameters.get());

    float stepU = 1.e-2;
    for (double u = 0; u < uMax + stepU; u += stepU) {
      for (double v = 0; v < uMax + stepU; v += stepU) {
        double f0 = F(u, v);
        float fSpline;
        spline.interpolate(1,(const float*)parameters.get(), u, v, &fSpline);
        statDf += (fSpline - f0) * (fSpline - f0);
        statN++;
      }
    }
    cout << "Spline2D standard deviation   : " << sqrt(statDf / statN) << std::endl;

    if (draw) {
      delete nt;
      delete knots;
      nt = new TNtuple("nt", "nt", "u:v:f0:fSpline");
      float stepU = 5.e-2;
      for (double u = 0; u < uMax + stepU; u += stepU) {
        for (double v = 0; v < uMax + stepU; v += stepU) {

          double f0 = F(u, v);
          float fSpline;
          spline.interpolate(1,(const float*)parameters.get(), u, v, &fSpline);
          nt->Fill(u, v, f0, fSpline);
        }
      }
      nt->SetMarkerStyle(8);

      nt->SetMarkerSize(.5);
      nt->SetMarkerColor(kBlue);
      nt->Draw("fSpline:u:v", "", "P");

      nt->SetMarkerColor(kGray);
      nt->SetMarkerSize(2.);
      nt->Draw("f0:u:v", "", "P,same");

      nt->SetMarkerSize(.5);
      nt->SetMarkerColor(kBlue);
      nt->Draw("fSpline:u:v", "", "P,same");

      knots = new TNtuple("knots", "knots", "type:u:v:f");
      for (int i = 0; i < nKnots; i++) {
        for (int j = 0; j < nKnots; j++) {
          double u = spline.getGridU().getKnot(i).u;
          double v = spline.getGridV().getKnot(j).u;
          float f;
          spline.interpolate(1,(const float*)parameters.get(), u, v, &f);
          knots->Fill(1, u, v, f);
          /*
          int naxU = 0, naxV = 0;
          double du = 0, dv = 0;
          if (i < nKnots - 1) {
            du = (spline.getGridU().getKnot(i + 1).u - u) / (nAxiliaryPoints + 1);
            naxU = nAxiliaryPoints;
          }
          if (j < nKnots - 1) {
            dv = (spline.getGridV().getKnot(j + 1).u - v) / (nAxiliaryPoints + 1);
            naxV = nAxiliaryPoints;
          }
          for (int ii = 0; ii <= naxU; ii++) {
            double uu = u + du * (ii);
            for (int jj = 0; jj <= naxV; jj++) {
              if (ii == 0 && jj == 0)
                continue;
              double vv = v + dv * (jj);
              float ff;
              spline.interpolate(1,(const float*)parameters.get(), uu, vv, &ff);
              knots->Fill(2, uu, ff);
            }
          }
          */
        }
      }

      knots->SetMarkerStyle(8);
      knots->SetMarkerSize(1.5);
      knots->SetMarkerColor(kRed);
      knots->SetMarkerSize(1.5);
      knots->Draw("f:u:v", "type==1", "same"); // compact
      /*
        knots->SetMarkerColor(kBlack);
        knots->SetMarkerSize(1.);
        knots->Draw("f:u:v", "type==2", "same"); // compact, axiliary points
        */
      if (!ask())
        break;
    }
  }
  //delete canv;
  //delete nt;
  //delete knots;

  statDf = sqrt(statDf / statN);
  cout << "\n std dev for Compact Spline   : " << statDf << std::endl;
  if (statDf < 0.1)
    cout << "Everything is fine" << endl;
  else {
    cout << "Something is wrong!!" << endl;
    return -2;
  }
  return 0;
}

#endif // GPUCA_GPUCODE
