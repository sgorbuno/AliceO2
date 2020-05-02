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
#include <iostream>
#include <chrono>

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

  const int Ndim = 3;

  const int Fdegree = 4;

  double Fcoeff[Ndim][4 * (Fdegree + 1) * (Fdegree + 1)];

  int nKnots = 16/*4*/;
  const int nAxiliaryPoints = 5;
  int uMax = nKnots * 3;

  auto F = [&](float u, float v, float Fuv[]) -> void {
    double uu = u * TMath::Pi() / uMax;
    double vv = v * TMath::Pi() / uMax;
    for (int dim = 0; dim < Ndim; dim++) {
      double f = 0; // Fcoeff[dim][0]/2;
      for (int i = 1; i <= Fdegree; i++) {
        double cosu = TMath::Cos(i * uu);
        double sinu = TMath::Sin(i * uu);
        for (int j = 1; j <= Fdegree; j++) {
          double *c = &(Fcoeff[dim][4 * (i * Fdegree + j)]);
          double cosv = TMath::Cos(j * vv);
          double sinv = TMath::Sin(j * vv);
          f += c[0] * cosu * cosv;
          f += c[1] * cosu * sinv;
          f += c[2] * sinu * cosv;
          f += c[3] * sinu * sinv;
        }
      }
      Fuv[dim] = f;
    }
  };

  TCanvas *canv = nullptr;
  TNtuple *nt = nullptr;
  TNtuple *knots = nullptr;

  auto ask = [&]() -> bool {
    if (!canv)
      return 0;
    canv->Update();
    cout << "type 'q ' to exit" << endl;
    std::string str;
    std::getline(std::cin, str);
    return (str != "q" && str != ".q");
  };

  std::cout << "Test 2D interpolation with the compact spline" << std::endl;

  int nTries = 10;

  if (draw) {
    canv = new TCanvas("cQA", "Spline1D  QA", 1200, 600);
    nTries = 10000;
  }

  long double statDf = 0;
  long double statN = 0;

  for (int seed = 1; seed < nTries + 1; seed++) {
    //cout << "next try.." << endl;

    gRandom->SetSeed(seed);

    for (int dim = 0; dim < Ndim; dim++) {
      for (int i = 0; i < 4 * (Fdegree + 1) * (Fdegree + 1); i++) {
        Fcoeff[dim][i] = gRandom->Uniform(-1, 1);
      }
    }

    o2::gpu::SplineHelper2D helper;
    o2::gpu::Spline2D spline;
    // spline.constructRegular(nKnots, nKnots);

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

      if (nKnots != spline.getGridU().getNumberOfKnots() ||
          nKnots != spline.getGridV().getNumberOfKnots()) {
        cout << "warning: n knots changed during the initialisation " << nKnots
             << " -> " << spline.getNumberOfKnots() << std::endl;
        continue;
      }
    } while (0);

    std::string err = FlatObject::stressTest(spline);
    if (!err.empty()) {
      cout << "error at FlatObject functionality: " << err << endl;
      return -1;
    } else {
      // cout << "flat object functionality is ok" << endl;
    }

    int err2 = helper.setSpline(spline, nAxiliaryPoints, nAxiliaryPoints);
    if (err2 != 0) {
      cout << "Error by spline construction: " << helper.getLastError()
           << std::endl;
      return -1;
    }
    std::unique_ptr<float[]> parameters(
        new float[spline.getNumberOfParameters(Ndim)]);

    std::unique_ptr<float[]> mapF(
        new float[helper.getNumberOfDataPoints() * Ndim]);

    for (int i = 0; i < spline.getNumberOfParameters(Ndim); i++)
      parameters[i] = 0.f;

    int nPointsU = helper.getNumberOfDataPointsU();
    int nPointsV = helper.getNumberOfDataPointsV();
    for (int ipv = 0; ipv < nPointsV; ipv++) {
      float v = helper.getHelperV().getDataPoint(ipv).u;
      for (int ipu = 0; ipu < nPointsU; ipu++) {
        float u = helper.getHelperU().getDataPoint(ipu).u;
        float Fuv[Ndim];
        F(u, v, Fuv);
        for (int dim = 0; dim < Ndim; dim++) {
          mapF[(ipv * nPointsU + ipu) * Ndim + dim] = Fuv[dim];
        }
      }
    }

    helper.constructParameters(Ndim, mapF.get(), parameters.get());

    double stepU = .1;
    for (double u = 0; u < uMax; u += stepU) {
      for (double v = 0; v < uMax; v += stepU) {
        float f[Ndim];
        F(u, v, f);
        float s[Ndim];
        spline.interpolate(Ndim, parameters.get(), u, v, s);
        for (int dim = 0; dim < Ndim; dim++) {
          statDf += (s[dim] - f[dim]) * (s[dim] - f[dim]);
        }
        statN += Ndim;
        // cout << u << " " << v << ": f " << f << " s " << s << " df "
        //   << s - f << " " << sqrt(statDf / statN) << std::endl;
      }
    }
    // cout << "Spline2D standard deviation   : " << sqrt(statDf / statN)
    //   << std::endl;

    if (draw) {
      delete nt;
      delete knots;
      nt = new TNtuple("nt", "nt", "u:v:f:s");
      knots = new TNtuple("knots", "knots", "type:u:v:s");

      double stepU = .1;
      for (double u = 0; u < uMax; u += stepU) {
        for (double v = 0; v < uMax; v += stepU) {
          float f[Ndim];
          F(u, v, f);
          float s[Ndim];
          spline.interpolate(Ndim, parameters.get(), u, v, s);
          nt->Fill(u, v, f[0], s[0]);
        }
      }
      nt->SetMarkerStyle(8);

      nt->SetMarkerSize(.5);
      nt->SetMarkerColor(kBlue);
      nt->Draw("s:u:v", "", "");

      nt->SetMarkerColor(kGray);
      nt->SetMarkerSize(2.);
      nt->Draw("f:u:v", "", "same");

      nt->SetMarkerSize(.5);
      nt->SetMarkerColor(kBlue);
      nt->Draw("s:u:v", "", "same");
      
      const Vc::short_v te;
      double u;
      double v;
      float s[Ndim];
      Vc::float_v sVec[Ndim];
      Vc::float_v uVec;
      Vc::float_v vVec;
      float sum = 0;
      int c;
      auto cache = Vc::int_v::IndexType::IndexesFromZero();

      /*
      std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
      for(size_t x = 0; x < 10000; x++){
      for (int i = 0; i < nKnots; i++) {
        for (int j = 0; j < nKnots; j++) {
          u = spline.getGridU().getKnot(i).u;
          v = spline.getGridV().getKnot(j).u;
          spline.interpolate(Ndim, parameters.get(), u, v, s);
          sum += s[0];
        }}}
      std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now(); 
      std::cout << sum << "\n";
      sum = 0;
      
      std::chrono::steady_clock::time_point beginVecHor = std::chrono::steady_clock::now();
      for(size_t x = 0; x < 10000; x++){
      for (int i = 0; i < nKnots; i++) {
        for (int j = 0; j < nKnots; j+=Vc::float_v::size()) {
          uVec = spline.getGridU().getKnot(Vc::int_v(i)).u;
          vVec = spline.getGridV().getKnot(cache + j).u;
          spline.interpolateVecHorizontal(Ndim, parameters.get(), uVec, vVec, sVec);
          for(size_t k = 0; k < Vc::float_v::size(); k++)
          { 
            sum += sVec[0][k];
          }
        }}}
      std::chrono::steady_clock::time_point endVecHor = std::chrono::steady_clock::now();     
      std::cout << sum << "\n";
      sum = 0;
      
      std::chrono::steady_clock::time_point beginVecHorV2 = std::chrono::steady_clock::now();
      for(size_t x = 0; x < 10000; x++){
      for (int i = 0; i < nKnots; i++) {
        for (int j = 0; j < nKnots; j+=Vc::float_v::size()) {
          uVec = spline.getGridU().getKnot(Vc::int_v(i)).u;
          vVec = spline.getGridV().getKnot(cache + j).u;
          spline.interpolateVecHorizontalV2(Ndim, parameters.get(), uVec, vVec, sVec);
          for(size_t k = 0; k < Vc::float_v::size(); k++)
          { 
            sum += sVec[0][k];
          }
        }}}
      std::chrono::steady_clock::time_point endVecHorV2 = std::chrono::steady_clock::now();     
      std::cout << sum << "\n";
      sum = 0;
      
      std::chrono::steady_clock::time_point beginVecHorV3 = std::chrono::steady_clock::now();
      for(size_t x = 0; x < 10000; x++){
      for (int i = 0; i < nKnots; i++) {
        for (int j = 0; j < nKnots; j+=Vc::float_v::size()) {
          uVec = spline.getGridU().getKnot(Vc::int_v(i)).u;
          vVec = spline.getGridV().getKnot(cache + j).u;
          spline.interpolateVecHorizontalV3(Ndim, parameters.get(), uVec, vVec, sVec);
          for(size_t k = 0; k < Vc::float_v::size(); k++)
          { 
            sum += sVec[0][k];
          }
        }}}
      std::chrono::steady_clock::time_point endVecHorV3 = std::chrono::steady_clock::now();     
      std::cout << sum << "\n";
      sum = 0;
      
      std::chrono::steady_clock::time_point beginVecAoV = std::chrono::steady_clock::now();
      for(size_t x = 0; x < 10000; x++){
      for (int i = 0; i < nKnots; i++) {
        for (int j = 0; j < nKnots; j++) {
          u = spline.getGridU().getKnot(i).u;
          v = spline.getGridV().getKnot(j).u;
          spline.interpolateVecVerticalAoV(Ndim, parameters.get(), u, v, s);
          sum += s[0];
        }}}
      std::chrono::steady_clock::time_point endVecAoV = std::chrono::steady_clock::now();     
      std::cout << sum << "\n";
      sum = 0;
      
      std::chrono::steady_clock::time_point beginVecSingleVector = std::chrono::steady_clock::now();
      for(size_t x = 0; x < 10000; x++){
      for (int i = 0; i < nKnots; i++) {
        for (int j = 0; j < nKnots; j++) {
          u = spline.getGridU().getKnot(i).u;
          v = spline.getGridV().getKnot(j).u;
          spline.interpolateVecVerticalSingleVector(Ndim, parameters.get(), u, v, s);
          sum += s[0];
        }}}
      std::chrono::steady_clock::time_point endVecSingleVector = std::chrono::steady_clock::now();
      std::cout << sum << "\n";
      sum = 0;
      */
      std::chrono::steady_clock::time_point beginVec = std::chrono::steady_clock::now();
      for(size_t x = 0; x < 10000; x++){
      for (int i = 0; i < nKnots; i++) {
        for (int j = 0; j < nKnots; j++) {
          u = spline.getGridU().getKnot(i).u;
          v = spline.getGridV().getKnot(j).u;
          spline.interpolateVecVertical(Ndim, parameters.get(), u, v, s);
          sum += s[0];
        }}}
      std::chrono::steady_clock::time_point endVec = std::chrono::steady_clock::now();
      std::cout << sum << "\n";
      sum = 0;
      /*
      std::chrono::steady_clock::time_point beginVecTest = std::chrono::steady_clock::now();
      for(size_t x = 0; x < 10000; x++){
      for (int i = 0; i < nKnots; i++) {
        for (int j = 0; j < nKnots; j++) {
          u = spline.getGridU().getKnot(i).u;
          v = spline.getGridV().getKnot(j).u;
          spline.interpolateVecVerticalTest(Ndim, parameters.get(), u, v, s);
          sum += s[0];
        }}}
      std::chrono::steady_clock::time_point endVecTest = std::chrono::steady_clock::now();
      std::cout << sum << "\n";
      sum = 0;

      std::chrono::steady_clock::time_point beginOpti = std::chrono::steady_clock::now();
      for(size_t x = 0; x < 10000; x++){
      for (int i = 0; i < nKnots; i++) {
        for (int j = 0; j < nKnots; j++) {
          u = spline.getGridU().getKnot(i).u;
          v = spline.getGridV().getKnot(j).u;
          spline.interpolateOptimizedScalar(Ndim, parameters.get(), u, v, s);
          sum += s[0];
        }}}
      std::chrono::steady_clock::time_point endOpti = std::chrono::steady_clock::now(); 
      std::cout << sum << "\n";
      */
      //std::cout << "Time Scalar = " << std::chrono::duration_cast<std::chrono::microseconds> (end - begin).count() << "[µs]" << std::endl;
      //std::cout << "Time Scalar-Opti = " << std::chrono::duration_cast<std::chrono::microseconds> (endOpti - beginOpti).count() << "[µs]" << std::endl;
      //std::cout << "Time Vec-Ver-AoV = " << std::chrono::duration_cast<std::chrono::microseconds> (endVecAoV - beginVecAoV).count() << "[µs]" << std::endl;
      //std::cout << "Time Vec-Ver-SingleVector = " << std::chrono::duration_cast<std::chrono::microseconds> (endVecSingleVector - beginVecSingleVector).count() << "[µs]" << std::endl;
      std::cout << "Time Vec-Ver-4Interpolations = " << std::chrono::duration_cast<std::chrono::microseconds> (endVec - beginVec).count() << "[µs]" << std::endl;
      //std::cout << "Time Vec-Ver-Test = " << std::chrono::duration_cast<std::chrono::microseconds> (endVecTest - beginVecTest).count() << "[µs]" << std::endl;
      //std::cout << "Time Vec-Hor = " << std::chrono::duration_cast<std::chrono::microseconds> (endVecHor - beginVecHor).count() << "[µs]" << std::endl;
      //std::cout << "Time Vec-HorV2 = " << std::chrono::duration_cast<std::chrono::microseconds> (endVecHorV2 - beginVecHorV2).count() << "[µs]" << std::endl;
      //std::cout << "Time Vec-HorV3 = " << std::chrono::duration_cast<std::chrono::microseconds> (endVecHorV3 - beginVecHorV3).count() << "[µs]" << std::endl;
    
     int counter = 0;
     Vc::float_v sVecHor[Ndim * 4];
     Vc::float_v sVecHorV3[Ndim * 4];
     for (int i = 0; i < nKnots; i++) {
       for (int j = 0; j < nKnots; j++) {
          if (counter ==  Vc::float_v::size())
            counter = 0;
          double u = spline.getGridU().getKnot(i).u;
          double v = spline.getGridV().getKnot(j).u;
          if(!counter)
          {
            uVec = spline.getGridU().getKnot(Vc::int_v(i)).u;
            vVec = spline.getGridU().getKnot(cache + j).u;
          }
          float s[Ndim];
          float sVecVertical[Ndim]; 
          float sVecVerticalTest[Ndim]; 
          spline.interpolate(Ndim, parameters.get(), u, v, s);
          spline.interpolateVecVertical(Ndim, parameters.get(), u, v, sVecVertical);
          spline.interpolateVecVerticalTest(Ndim, parameters.get(), u, v, sVecVerticalTest);
          if(!counter)
          {
            spline.interpolateVecHorizontal(Ndim, parameters.get(), uVec, vVec, sVecHor);
            spline.interpolateVecHorizontalV3(Ndim, parameters.get(), uVec, vVec, sVecHorV3);
          }
          std::cout << "Interpolate:                  u = " << u << " v = " << v << " s = " << s[0] << "\n";
          std::cout << "InterpolateVecVertical:       u = " << u << " v = " << v << " s = " << sVecVertical[0] << "\n";
          std::cout << "InterpolateVecVerticalTest:   u = " << u << " v = " << v << " s = " << sVecVerticalTest[0] << "\n";
          std::cout << "InterpolateVecHorizontal:     u = " << uVec[counter] << " v = " << vVec[counter] << " s = " << sVecHor[0][counter] << "\n";
          std::cout << "InterpolateVecHorizontalV3:   u = " << uVec[counter] << " v = " << vVec[counter] << " s = " << sVecHorV3[0][counter] << "\n";

          counter++;
          knots->Fill(1, u, v, s[0]);
          /*
          int naxU = 0, naxV = 0;
          double du = 0, dv = 0;
          if (i < nKnots - 1) {
            du = (spline.getGridU().getKnot(i + 1).u - u) / (nAxiliaryPoints +
          1); naxU = nAxiliaryPoints;
          }
          if (j < nKnots - 1) {
            dv = (spline.getGridV().getKnot(j + 1).u - v) / (nAxiliaryPoints +
          1); naxV = nAxiliaryPoints;
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
      knots->Draw("s:u:v", "type==1", "same"); // compact
      /*
        knots->SetMarkerColor(kBlack);
        knots->SetMarkerSize(1.);
        knots->Draw("f:u:v", "type==2", "same"); // compact, axiliary points
        */
      if (!ask())
        break;
    }
  }
  // delete canv;
  // delete nt;
  // delete knots;

  statDf = sqrt(statDf / statN);
  cout << "\n std dev for Compact Spline   : " << statDf << std::endl;
  if (statDf < 0.15)
    cout << "Everything is fine" << endl;
  else {
    cout << "Something is wrong!!" << endl;
    return -2;
  }
  return 0;
}

#endif // GPUCA_GPUCODE
