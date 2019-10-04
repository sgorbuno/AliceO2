// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  CompactSplineIrregular1D.cxx
/// \brief Implementation of CompactSplineIrregular1D class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#include "CompactSplineIrregular1D.h"
#include <cmath>

#if !defined(GPUCA_GPUCODE) // code invisible on GPU
#include <vector>
#include <algorithm>
#include <iostream>
#endif

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE) // code invisible on GPU and in the standalone compilation
#include "TRandom.h"
#include "Riostream.h"
#include "TMath.h"
#include "CompactSplineHelper.h"
#include "TCanvas.h"
#include "TNtuple.h"
#endif

using namespace GPUCA_NAMESPACE::gpu;

CompactSplineIrregular1D::CompactSplineIrregular1D() : FlatObject(), mNumberOfKnots(0), mUmax(0), mBin2KnotMap(0)
{
  /// Default constructor. Creates an empty uninitialised object
}

void CompactSplineIrregular1D::destroy()
{
  /// See FlatObject for description
  mNumberOfKnots = 0;
  mUmax = 0;
  mBin2KnotMap = nullptr;
  FlatObject::destroy();
}

#if !defined(GPUCA_GPUCODE)
void CompactSplineIrregular1D::cloneFromObject(const CompactSplineIrregular1D& obj, char* newFlatBufferPtr)
{
  /// See FlatObject for description
  const char* oldFlatBufferPtr = obj.mFlatBufferPtr;
  FlatObject::cloneFromObject(obj, newFlatBufferPtr);
  mNumberOfKnots = obj.mNumberOfKnots;
  mUmax = obj.mUmax;
  mBin2KnotMap = FlatObject::relocatePointer(oldFlatBufferPtr, mFlatBufferPtr, obj.mBin2KnotMap);
}

void CompactSplineIrregular1D::moveBufferTo(char* newFlatBufferPtr)
{
  /// See FlatObject for description
  const char* oldFlatBufferPtr = mFlatBufferPtr;
  FlatObject::moveBufferTo(newFlatBufferPtr);
  mBin2KnotMap = FlatObject::relocatePointer(oldFlatBufferPtr, mFlatBufferPtr, mBin2KnotMap);
}
#endif

void CompactSplineIrregular1D::setActualBufferAddress(char* actualFlatBufferPtr)
{
  /// See FlatObject for description
  mBin2KnotMap = FlatObject::relocatePointer(mFlatBufferPtr, actualFlatBufferPtr, mBin2KnotMap);
  FlatObject::setActualBufferAddress(actualFlatBufferPtr);
}

void CompactSplineIrregular1D::setFutureBufferAddress(char* futureFlatBufferPtr)
{
  /// See FlatObject for description
  mBin2KnotMap = FlatObject::relocatePointer(mFlatBufferPtr, futureFlatBufferPtr, mBin2KnotMap);
  FlatObject::setFutureBufferAddress(futureFlatBufferPtr);
}

#if !defined(GPUCA_GPUCODE)

void CompactSplineIrregular1D::construct(int numberOfKnots, const int inputKnots[])
{
  /// Constructor
  ///
  /// Number of created knots may differ from the input values:
  /// - Edge knots {0} and {Umax} will be added if they are not present.
  /// - Duplicated knots, knots with a negative coordinate will be deleted
  /// - At least 2 knots will be created
  ///
  /// \param numberOfKnots     Number of knots in knots[] array
  /// \param knots             Array of knot positions (integer values)
  ///

  FlatObject::startConstruction();

  std::vector<int> knotU;

  { // reorganize knots

    std::vector<int> tmp;
    for (int i = 0; i < numberOfKnots; i++)
      tmp.push_back(inputKnots[i]);
    std::sort(tmp.begin(), tmp.end());

    knotU.push_back(0); // obligatory knot at 0.0

    for (int i = 0; i < numberOfKnots; ++i) {
      if (knotU.back() < tmp[i])
        knotU.push_back(tmp[i]);
    }
    if (knotU.back() < 1)
      knotU.push_back(1);
  }

  mNumberOfKnots = knotU.size();
  mUmax = knotU.back();
  int bin2KnotMapOffset = mNumberOfKnots * sizeof(CompactSplineIrregular1D::Knot);

  FlatObject::finishConstruction(bin2KnotMapOffset + (mUmax + 1) * sizeof(int));

  mBin2KnotMap = reinterpret_cast<int*>(mFlatBufferPtr + bin2KnotMapOffset);

  CompactSplineIrregular1D::Knot* s = getKnotsNonConst();

  for (int i = 0; i < mNumberOfKnots; i++) {
    s[i].u = knotU[i];
  }

  for (int i = 0; i < mNumberOfKnots - 1; i++) {
    s[i].Li = 1. / (s[i + 1].u - s[i].u); // do division in double
  }

  s[mNumberOfKnots - 1].Li = 0.f; // the value will not be used, we define it for consistency

  // Set up map (U bin) -> (knot index)

  int* map = getBin2KnotMapNonConst();

  int iKnotMax = mNumberOfKnots - 2;

  //
  // With iKnotMax=nKnots-2 we map the U==Umax coordinate to the [nKnots-2, nKnots-1] segment.
  // This trick allows one to avoid a special condition for this edge case.
  // Any U from [0,Umax] is mapped to some knot i such, that the knot i+1 is always exist
  //

  for (int u = 0, iKnot = 0; u <= mUmax; u++) {
    if ((knotU[iKnot + 1] == u) && (iKnot < iKnotMax)) {
      iKnot = iKnot + 1;
    }
    map[u] = iKnot;
  }
}

void CompactSplineIrregular1D::constructRegular(int numberOfKnots)
{
  /// Constructor for a regular spline
  /// \param numberOfKnots     Number of knots
  ///

  if (numberOfKnots < 2)
    numberOfKnots = 2;

  std::vector<int> knots(numberOfKnots);
  for (int i = 0; i < numberOfKnots; i++) {
    knots[i] = i;
  }
  construct(numberOfKnots, knots.data());
}
#endif

void CompactSplineIrregular1D::print() const
{
#if !defined(GPUCA_GPUCODE)
  std::cout << " Compact Spline 1D: " << std::endl;
  std::cout << "  mNumberOfKnots = " << mNumberOfKnots << std::endl;
  std::cout << "  mUmax = " << mUmax << std::endl;
  std::cout << "  mBin2KnotMap = " << (void*)mBin2KnotMap << std::endl;
  std::cout << "  knots: ";
  for (int i = 0; i < mNumberOfKnots; i++) {
    std::cout << getKnot(i).u << " ";
  }
  std::cout << std::endl;
#endif
}

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE) // code invisible on GPU and in the standalone compilation

int CompactSplineIrregular1D::test(bool draw)
{
  using namespace std;

  const int funcN = 10;
  double funcC[2 * funcN + 2];

  int nKnots = 4;
  const int nAxiliaryPoints = 5;
  int uMax = nKnots * 3;

  auto F = [&](float u) -> float {
    double uu = u * TMath::Pi() / uMax;
    double f = 0; //funcC[0]/2;
    for (int i = 1; i <= funcN; i++) {
      f += funcC[2 * i] * TMath::Cos(i * uu) + funcC[2 * i + 1] * TMath::Sin(i * uu);
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
    << "Test 1D interpolation with the compact spline" << std::endl;

  int nTries = 100;

  if (draw) {
    canv = new TCanvas("cQA", "CompactSplineIrregular1D  QA", 2000, 1000);
    nTries = 10000;
  }

  double statDf = 0;
  double statN = 0;

  for (int seed = 1; seed < nTries; seed++) {

    gRandom->SetSeed(seed);

    for (int i = 0; i <= funcN; i++) {
      funcC[i] = gRandom->Uniform(-1, 1);
    }

    CompactSplineHelper helper;
    CompactSplineIrregular1D spline;

    int knotsU[nKnots];
    do {
      knotsU[0] = 0;
      double du = 1. * uMax / (nKnots - 1);
      for (int i = 1; i < nKnots; i++) {
        knotsU[i] = (int)(i * du); //+ gRandom->Uniform(-du / 3, du / 3);
      }
      knotsU[nKnots - 1] = uMax;
      spline.construct(nKnots, knotsU);

      if (nKnots != spline.getNumberOfKnots()) {
        cout << "warning: n knots changed during the initialisation " << nKnots << " -> " << spline.getNumberOfKnots() << std::endl;
        continue;
      }
    } while (0);

    std::string err = FlatObject::stressTest(spline);
    if (!err.empty()) {
      cout << "error at FlatObject functionality: " << err << endl;
      return -1;
    } else {
      cout << "flat object functionality is ok" << endl;
    }

    nKnots = spline.getNumberOfKnots();
    cout << "mark 0: nKnots = " << nKnots << endl;
    std::unique_ptr<float[]> data = helper.create(spline, F, nAxiliaryPoints);
    cout << "mark 1" << endl;
    if (data == nullptr) {
      cout << "can not create data array for the spline" << endl;
      return -3;
    }

    float stepU = 1.e-2;
    for (double u = 0; u < uMax + stepU; u += stepU) {
      double f0 = F(u);
      double fSpline = spline.getSpline((const float*)data.get(), u);
      statDf += (fSpline - f0) * (fSpline - f0);
      statN++;
    }
    //cout << "std dev Compact   : " << sqrt(statDf / statN) << std::endl;

    if (draw) {
      delete nt;
      delete knots;
      nt = new TNtuple("nt", "nt", "u:f0:fSpline");
      float stepU = 1.e-4;
      for (double u = 0; u < uMax + stepU; u += stepU) {
        double f0 = F(u);
        double fSpline = spline.getSpline((const float*)data.get(), u);
        nt->Fill(u, f0, fSpline);
      }

      nt->SetMarkerStyle(8);

      nt->SetMarkerSize(.5);
      nt->SetMarkerColor(kBlue);
      nt->Draw("fSpline:u", "", "P");

      nt->SetMarkerColor(kGray);
      nt->SetMarkerSize(2.);
      nt->Draw("f0:u", "", "P,same");

      nt->SetMarkerSize(.5);
      nt->SetMarkerColor(kBlue);
      nt->Draw("fSpline:u", "", "P,same");

      knots = new TNtuple("knots", "knots", "type:u:f");
      for (int i = 0; i < nKnots; i++) {
        double u = spline.getKnot(i).u;
        double f = spline.getSpline((const float*)data.get(), u);
        knots->Fill(1, u, f);
        if (i < nKnots - 1) {
          double u1 = spline.getKnot(i + 1).u;
          int nax = nAxiliaryPoints;
          double du = (u1 - u) / (nax + 1);
          for (int j = 0; j < nax; j++) {
            double uu = u + du * (j + 1);
            double ff = spline.getSpline((const float*)data.get(), uu);
            knots->Fill(2, uu, ff);
          }
        }
      }

      knots->SetMarkerStyle(8);
      knots->SetMarkerSize(1.5);
      knots->SetMarkerColor(kRed);
      knots->SetMarkerSize(1.5);
      knots->Draw("f:u", "type==1", "same"); // compact
      knots->SetMarkerColor(kBlack);
      knots->SetMarkerSize(1.);
      knots->Draw("f:u", "type==2", "same"); // compact, axiliary points

      if (!ask())
        break;
    }
  }
  delete canv;
  delete nt;
  delete knots;

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