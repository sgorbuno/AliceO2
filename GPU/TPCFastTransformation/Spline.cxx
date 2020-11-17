// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  Spline.cxx
/// \brief Implementation of Spline class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#include "Spline.h"

#if !defined(GPUCA_GPUCODE)
#include <iostream>
#endif

#if !defined(GPUCA_ALIGPUCODE) && !defined(GPUCA_STANDALONE) // code invisible on GPU and in the standalone compilation
#include "TRandom.h"
#include "Riostream.h"
#include "TMath.h"
#include "SplineHelper.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "TFile.h"

templateClassImp(GPUCA_NAMESPACE::gpu::SplineBase);

#endif

#include "Spline2D.h"

using namespace std;
using namespace GPUCA_NAMESPACE::gpu;

template <typename DataT, bool isConsistentT>
SplineBase<DataT, isConsistentT>::SplineBase(int nXdim, int nFdim)
  : FlatObject(), mXdim(nXdim), mFdim(nFdim), mNparameters(0), mGrid(nullptr), mFparameters(nullptr)
{
  recreate(nullptr, nullptr);
}

template <typename DataT, bool isConsistentT>
void SplineBase<DataT, isConsistentT>::destroy()
{
  /// See FlatObject for description
  mXdim = 0;
  mFdim = 0;
  mNparameters = 0;
  mGrid = nullptr;
  mFparameters = nullptr;
  FlatObject::destroy();
}

template <typename DataT, bool isConsistentT>
void SplineBase<DataT, isConsistentT>::setActualBufferAddress(char* actualFlatBufferPtr)
{
  /// See FlatObject for description

  FlatObject::setActualBufferAddress(actualFlatBufferPtr);

  mGrid = reinterpret_cast<Spline1D<DataT>*>(mFlatBufferPtr);

  int offset = sizeof(*mGrid) * mXdim;

  for (int i = 0; i < mXdim; i++) {
    offset = alignSize(offset, mGrid[i].getBufferAlignmentBytes());
    mGrid[i].setActualBufferAddress(mFlatBufferPtr + offset);
    offset += mGrid[i].getFlatBufferSize();
  }

  if (isConsistentT) {
    offset = alignSize(offset, getParameterAlignmentBytes());
    mFparameters = reinterpret_cast<DataT*>(mFlatBufferPtr + offset);
    //offset += getSizeOfParameters();
  } else {
    mFparameters = nullptr;
  }
}

template <typename DataT, bool isConsistentT>
void SplineBase<DataT, isConsistentT>::setFutureBufferAddress(char* futureFlatBufferPtr)
{
  /// See FlatObject for description

  if (isConsistentT) {
    mFparameters = relocatePointer(mFlatBufferPtr, futureFlatBufferPtr, mFparameters);
  } else {
    mFparameters = nullptr;
  }
  for (int i = 0; i < mXdim; i++) {
    char* buffer = relocatePointer(mFlatBufferPtr, futureFlatBufferPtr, mGrid[i].getFlatBufferPtr());
    mGrid[i].setFutureBufferAddress(buffer);
  }
  mGrid = relocatePointer(mFlatBufferPtr, futureFlatBufferPtr, mGrid);
  FlatObject::setFutureBufferAddress(futureFlatBufferPtr);
}

template <typename DataT, bool isConsistentT>
void SplineBase<DataT, isConsistentT>::print() const
{
  printf(" Irregular Spline %dD->%dD: \n", mXdim, mFdim);
  for (int i = 0; i < mXdim; i++) {
    printf(" grid U%d: \n", i);
    mGrid[i].print();
  }
}

#if !defined(GPUCA_GPUCODE)
template <typename DataT, bool isConsistentT>
void SplineBase<DataT, isConsistentT>::cloneFromObject(const SplineBase<DataT, isConsistentT>& obj, char* newFlatBufferPtr)
{
  /// See FlatObject for description
  if (mXdim != obj.mXdim || isConsistentT && mFdim != obj.mFdim) {
    assert(0);
    return;
  }

  const char* oldFlatBufferPtr = obj.mFlatBufferPtr;

  FlatObject::cloneFromObject(obj, newFlatBufferPtr);
  mGrid = FlatObject::relocatePointer(oldFlatBufferPtr, mFlatBufferPtr, obj.mGrid);
  for (int i = 0; i < mXdim; i++) {
    char* buffer = FlatObject::relocatePointer(oldFlatBufferPtr, mFlatBufferPtr, obj.mGrid[i].getFlatBufferPtr());
    mGrid[i].cloneFromObject(obj.mGrid[i], buffer);
  }

  if (isConsistentT) {
    mFparameters = FlatObject::relocatePointer(oldFlatBufferPtr, mFlatBufferPtr, obj.mFparameters);
  } else {
    mFparameters = nullptr;
  }
}

template <typename DataT, bool isConsistentT>
void SplineBase<DataT, isConsistentT>::moveBufferTo(char* newFlatBufferPtr)
{
  /// See FlatObject for description
  char* oldFlatBufferPtr = mFlatBufferPtr;
  FlatObject::moveBufferTo(newFlatBufferPtr);
  char* currFlatBufferPtr = mFlatBufferPtr;
  mFlatBufferPtr = oldFlatBufferPtr;
  setActualBufferAddress(currFlatBufferPtr);
}

template <typename DataT, bool isConsistentT>
void SplineBase<DataT, isConsistentT>::recreate(
  const int numberOfKnots[/* mXdim */], const int* const knots[/* mXdim */])
{
  /// Constructor for an irregular spline

  FlatObject::startConstruction();

  Spline1D<DataT> vGrids[mXdim];

  mNparameters = getNumberOfParametersPerKnot();
  for (int i = 0; i < mXdim; i++) {
    if (knots) {
      vGrids[i].recreate(numberOfKnots[i], knots[i], 0);
    } else if (numberOfKnots) {
      vGrids[i].recreate(numberOfKnots[i], 0);
    } else {
      vGrids[i].recreate(2, 0);
    }
    mNparameters *= vGrids[i].getNumberOfKnots();
  }

  int offset = sizeof(Spline1D<DataT>) * mXdim;

  for (int i = 0; i < mXdim; i++) {
    offset = alignSize(offset, vGrids[i].getBufferAlignmentBytes());
    offset += vGrids[i].getFlatBufferSize();
  }

  if (isConsistentT) {
    offset = alignSize(offset, getParameterAlignmentBytes());
    offset += getSizeOfParameters();
  }

  FlatObject::finishConstruction(offset);

  mGrid = reinterpret_cast<Spline1D<DataT>*>(mFlatBufferPtr);

  offset = sizeof(Spline1D<DataT>) * mXdim;

  for (int i = 0; i < mXdim; i++) {
    new (&mGrid[i]) Spline1D<DataT>; // constructor
    offset = alignSize(offset, mGrid[i].getBufferAlignmentBytes());
    mGrid[i].cloneFromObject(vGrids[i], mFlatBufferPtr + offset);
    offset += mGrid[i].getFlatBufferSize();
  }

  if (isConsistentT) {
    offset = alignSize(offset, getParameterAlignmentBytes());
    mFparameters = reinterpret_cast<DataT*>(mFlatBufferPtr + offset);
    offset += getSizeOfParameters();

    for (int i = 0; i < getNumberOfParameters(); i++) {
      mFparameters[i] = 0;
    }
  } else {
    mFparameters = nullptr;
  }
}

template <typename DataT, bool isConsistentT>
void SplineBase<DataT, isConsistentT>::recreate(
  const int numberOfKnots[/* mXdim */])
{
  /// Constructor for a regular spline
  recreate(numberOfKnots, nullptr);
}
#endif

#if !defined(GPUCA_ALIGPUCODE) && !defined(GPUCA_STANDALONE) // code invisible on GPU and in the standalone compilation

template <typename DataT, bool isConsistentT>
void SplineBase<DataT, isConsistentT>::
  approximateFunction(
    const DataT xMin[/* mXdim */], const DataT xMax[/* mXdim */],
    std::function<void(const DataT x[/* mXdim */], DataT f[/* mFdim */])> F,
    const int nAxiliaryDataPoints[/* mXdim */])
{
  /// approximate a function F with this spline
  SplineHelper<DataT> helper;
  helper.approximateFunction(*this, xMin, xMax, F, nAxiliaryDataPoints);
}

template <typename DataT, bool isConsistentT>
int SplineBase<DataT, isConsistentT>::writeToFile(TFile& outf, const char* name)
{
  /// write a class object to the file
  return FlatObject::writeToFile(*this, outf, name);
}

template <typename DataT, bool isConsistentT>
SplineBase<DataT, isConsistentT>* SplineBase<DataT, isConsistentT>::readFromFile(
  TFile& inpf, const char* name)
{
  /// read a class object from the file
  return FlatObject::readFromFile<SplineBase<DataT, isConsistentT>>(inpf, name);
}

// TODO: Implement
/*
  using namespace std;

  const int Ndim = 3;

  const int Fdegree = 4;

  double Fcoeff[Ndim][4 * (Fdegree + 1) * (Fdegree + 1)];

  int nKnots = 4;
  const int nAxiliaryPoints = 1;
  int uMax = nKnots * 3;

  auto F = [&](DataT u, DataT v, DataT Fuv[]) {
    double uu = u * TMath::Pi() / uMax;
    double vv = v * TMath::Pi() / uMax;
    for (int dim = 0; dim < Ndim; dim++) {
      double f = 0; // Fcoeff[dim][0]/2;
      for (int i = 1; i <= Fdegree; i++) {
        double cosu = TMath::Cos(i * uu);
        double sinu = TMath::Sin(i * uu);
        for (int j = 1; j <= Fdegree; j++) {
          double* c = &(Fcoeff[dim][4 * (i * Fdegree + j)]);
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

  TCanvas* canv = nullptr;
  TNtuple* nt = nullptr;
  TNtuple* knots = nullptr;

  auto ask = [&]() -> bool {
    if (!canv) {
      return 0;
    }
    canv->Update();
    cout << "type 'q ' to exit" << endl;
    std::string str;
    std::getline(std::cin, str);
    return (str != "q" && str != ".q");
  };

  std::cout << "Test  interpolation with the compact spline" << std::endl;

  int nTries = 10;

  if (draw) {
    canv = new TCanvas("cQA", "Spline  QA", 1500, 800);
    nTries = 10000;
  }

  long double statDf = 0;
  long double statDf1D = 0;
  long double statN = 0;

  for (int seed = 1; seed < nTries + 1; seed++) {
    //cout << "next try.." << endl;

    gRandom->SetSeed(seed);

    for (int dim = 0; dim < Ndim; dim++) {
      for (int i = 0; i < 4 * (Fdegree + 1) * (Fdegree + 1); i++) {
        Fcoeff[dim][i] = gRandom->Uniform(-1, 1);
      }
    }

    Spline<DataT, Ndim> spline;

    int knotsU[nKnots], knotsV[nKnots];
    do {
      knotsU[0] = 0;
      knotsV[0] = 0;
      double du = 1. * uMax / (nKnots - 1);
      for (int i = 1; i < nKnots; i++) {
        knotsU[i] = (int)(i * du); // + gRandom->Uniform(-du / 3, du / 3);
        knotsV[i] = (int)(i * du); // + gRandom->Uniform(-du / 3, du / 3);
      }
      knotsU[nKnots - 1] = uMax;
      knotsV[nKnots - 1] = uMax;
      spline.recreate(nKnots, knotsU, nKnots, knotsV);

      if (nKnots != spline.getGridU1().getNumberOfKnots() ||
          nKnots != spline.getGridU2().getNumberOfKnots()) {
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

    // Ndim-D spline
    spline.approximateFunction(0., uMax, 0., uMax, F, 4, 4);

    //if (itry == 0)
    if (1) {
      TFile outf("testSpline.root", "recreate");
      if (outf.IsZombie()) {
        cout << "Failed to open output file testSpline.root " << std::endl;
      } else {
        const char* name = "splinetest";
        spline.writeToFile(outf, name);
        Spline<DataT, Ndim>* p = Spline<DataT, Ndim>::readFromFile(outf, name);
        if (p == nullptr) {
          cout << "Failed to read Spline1D from file testSpline1D.root " << std::endl;
        } else {
          spline = *p;
        }
        outf.Close();
      }
    }

    // 1-D splines for each of Ndim dimensions

    Spline<DataT, 1> splines1D[Ndim];

    for (int dim = 0; dim < Ndim; dim++) {
      auto F1 = [&](DataT x1, DataT x2, DataT f[]) {
        DataT ff[Ndim];
        F(x1, x2, ff);
        f[0] = ff[dim];
      };
      splines1D[dim].recreate(nKnots, knotsU, nKnots, knotsV);
      splines1D[dim].approximateFunction(0., uMax, 0., uMax, F1, 4, 4);
    }

    double stepU = .1;
    for (double u = 0; u < uMax; u += stepU) {
      for (double v = 0; v < uMax; v += stepU) {
        DataT f[Ndim];
        F(u, v, f);
        DataT s[Ndim];
        DataT s1;
        spline.interpolate(u, v, s);
        for (int dim = 0; dim < Ndim; dim++) {
          statDf += (s[dim] - f[dim]) * (s[dim] - f[dim]);
          splines1D[dim].interpolate(u, v, &s1);
          statDf1D += (s[dim] - s1) * (s[dim] - s1);
        }
        statN += Ndim;
        // cout << u << " " << v << ": f " << f << " s " << s << " df "
        //   << s - f << " " << sqrt(statDf / statN) << std::endl;
      }
    }
    // cout << "Spline standard deviation   : " << sqrt(statDf / statN)
    //   << std::endl;

    if (draw) {
      delete nt;
      delete knots;
      nt = new TNtuple("nt", "nt", "u:v:f:s");
      knots = new TNtuple("knots", "knots", "type:u:v:s");
      double stepU = .3;
      for (double u = 0; u < uMax; u += stepU) {
        for (double v = 0; v < uMax; v += stepU) {
          DataT f[Ndim];
          F(u, v, f);
          DataT s[Ndim];
          spline.interpolate(u, v, s);
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

      for (int i = 0; i < nKnots; i++) {
        for (int j = 0; j < nKnots; j++) {
          double u = spline.getGridU1().getKnot(i).u;
          double v = spline.getGridU2().getKnot(j).u;
          DataT s[Ndim];
          spline.interpolate(u, v, s);
          knots->Fill(1, u, v, s[0]);
        }
      }

      knots->SetMarkerStyle(8);
      knots->SetMarkerSize(1.5);
      knots->SetMarkerColor(kRed);
      knots->SetMarkerSize(1.5);
      knots->Draw("s:u:v", "type==1", "same"); // knots

      if (drawDataPoints) {
        SplineHelper<DataT> helper;
        helper.setSpline(spline, 4, 4);
        for (int ipu = 0; ipu < helper.getHelperU1().getNumberOfDataPoints(); ipu++) {
          const typename SplineHelper1D<DataT>::DataPoint& pu = helper.getHelperU1().getDataPoint(ipu);
          for (int ipv = 0; ipv < helper.getHelperU2().getNumberOfDataPoints(); ipv++) {
            const typename SplineHelper1D<DataT>::DataPoint& pv = helper.getHelperU2().getDataPoint(ipv);
            if (pu.isKnot && pv.isKnot) {
              continue;
            }
            DataT s[Ndim];
            spline.interpolate(pu.u, pv.u, s);
            knots->Fill(2, pu.u, pv.u, s[0]);
          }
        }
        knots->SetMarkerColor(kBlack);
        knots->SetMarkerSize(1.);
        knots->Draw("s:u:v", "type==2", "same"); // data points
      }

      if (!ask()) {
        break;
      }
    }
  }
  // delete canv;
  // delete nt;
  // delete knots;

  statDf = sqrt(statDf / statN);
  statDf1D = sqrt(statDf1D / statN);

  cout << "\n std dev for Spline   : " << statDf << std::endl;
  cout << " mean difference between 1-D and " << Ndim
       << "-D splines   : " << statDf1D << std::endl;

  if (statDf < 0.15 && statDf1D < 1.e-20) {
    cout << "Everything is fine" << endl;
  } else {
    cout << "Something is wrong!!" << endl;
    return -2;
  }
*/

#ifdef XXX
////////////////
////TESTFUNCTION 2D
template <typename DataT, bool isConsistentT>
int SplineBase<DataT, isConsistentT>::test(const bool draw, const bool drawDataPoints)
{
  using namespace std;

  const int Ndim = 3;
  const int Fdegree = 4;
  double Fcoeff[Ndim][4 * (Fdegree + 1) * (Fdegree + 1)];

  constexpr int nKnots = 4;
  constexpr int nAxiliaryPoints = 1;
  constexpr int uMax = nKnots * 3;
  auto F = [&](DataT u[], DataT Fuv[]) {
    const double scale = TMath::Pi() / uMax;
    double uu = u[0] * scale;
    double vv = u[1] * scale;
    double cosu[Fdegree + 1], sinu[Fdegree + 1], cosv[Fdegree + 1], sinv[Fdegree + 1];
    double ui = 0, vi = 0;
    for (int i = 0; i <= Fdegree; i++, ui += uu, vi += vv) {
      cosu[i] = cos(ui);
      sinu[i] = sin(ui);
      cosv[i] = cos(vi);
      sinv[i] = sin(vi);
    }
    for (int dim = 0; dim < Ndim; dim++) {
      double f = 0; // Fcoeff[dim][0]/2;
      for (int i = 1; i <= Fdegree; i++) {
        for (int j = 1; j <= Fdegree; j++) {
          double* c = &(Fcoeff[dim][4 * (i * Fdegree + j)]);
          f += c[0] * cosu[i] * cosv[j];
          f += c[1] * cosu[i] * sinv[j];
          f += c[2] * sinu[i] * cosv[j];
          f += c[3] * sinu[i] * sinv[j];
        }
      }
      Fuv[dim] = f;
    }
  };

  int seed = 1;
  gRandom->SetSeed(seed);

  for (int dim = 0; dim < Ndim; dim++) {
    for (int i = 0; i < 4 * (Fdegree + 1) * (Fdegree + 1); i++) {
      Fcoeff[dim][i] = gRandom->Uniform(-1, 1);
    }
  }
  std::cout << "Fcoeff: " << std::endl;
  for (int dim = 0; dim < Ndim; dim++) {
    for (int i = 0; i < 4 * (Fdegree + 1) * (Fdegree + 1); i++) {
      std::cout << Fcoeff[dim][i] << ", " << std::endl;
    }
  }
  std::cout << std::endl;

  TCanvas* canv = nullptr;
  TNtuple* nt = nullptr;
  TNtuple* knots = nullptr;
  /* 
  auto ask = [&]() -> bool {
    if (!canv) {
      return 0;
    }
    canv->Update();
    cout << "type 'q ' to exit" << endl;
    std::string str;
    std::getline(std::cin, str);
    return (str != "q" && str != ".q");
  };
  std::cout << "Test 2D interpolation with the compact  ND spline" << std::endl;
  int nTries = 10;
  if (draw) {
    canv = new TCanvas("cQA", "Spline2D  QA", 1500, 800);
    nTries = 10000;
  }

  long double statDf = 0;
  long double statDf1D = 0;
  long double statN = 0;

  for (int seed = 1; seed < nTries + 1; seed++) {
    //cout << "next try.." << endl;

    gRandom->SetSeed(seed);

    for (int dim = 0; dim < Ndim; dim++) {
      for (int i = 0; i < 4 * (Fdegree + 1) * (Fdegree + 1); i++) {
        Fcoeff[dim][i] = gRandom->Uniform(-1, 1);
      }
    }
    //EINGEFÜGT: weil wegen instanziierung anders
    const int NdimX = 2;
    const int nDimF = 1;
    
    int **knotsU = new int* [nDimX];
    
    do {
      knotsU[0][0] = 0;
      knotsU[1][0] = 0;
      double du = 1. * uMax / (nKnots - 1);
      for (int i = 1; i < nKnots; i++) {
        knotsU[0][i] = (int)(i * du); // + gRandom->Uniform(-du / 3, du / 3);
        knotsU[1][i] = (int)(i * du); // + gRandom->Uniform(-du / 3, du / 3);
      }
      knotsU[0][nKnots - 1] = uMax;
      knotsU[1][nKnots - 1] = uMax;
      Spline<double, NdimX, nDimF> spline(nKnots, knotsU);
      //WARNUNG AUSGELASSENax;
      //spline.recreate(nKnots, kn
    } while (0);
    // FLAT OBJECT STRESSTEST AUSGELASSEN
    double xMin[NdimX]={0.,0.};
    double xMax[NdimX]= {3.,3.};
    int nAxiliaryDataPoints[NdimX]={4,4};
    spline.approximateFunction(xMin, xMax, F, nAxiliaryDataPoints);
    
    //WRITE TO TESTFILE AUSGELASSEN

    //1D SPLINES AUSGELASSEN
    //Standard derivation ausgelassen


    if (draw) {
      delete nt;
      delete knots;
      nt = new TNtuple("nt", "nt", "u:v:f:s");
      knots = new TNtuple("knots", "knots", "type:u:v:s");
      double stepU = .3;
      for (double u = 0; u < uMax; u += stepU) {
        for (double v = 0; v < uMax; v += stepU) {
          DataT f[Ndim];
          DataT inputx[NdimX];
          inputx[0] =u;
          inputx[1] = v;
          F(inputx, f);
          DataT s[Ndim];
          spline.interpolate(inputx, s);
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

      for (int i = 0; i < nKnots; i++) {
        for (int j = 0; j < nKnots; j++) {
          double inputx[2];
          inputx[0] = spline.getGrid(0).getKnot(i).u;
          inputx[1] = spline.getGrid(1).getKnot(j).u;
          DataT s[Ndim];
          spline.interpolate(inputx, s);
          knots->Fill(1, inputx[0], inputx[1], s[0]);
        }
      }

      knots->SetMarkerStyle(8);
      knots->SetMarkerSize(1.5);
      knots->SetMarkerColor(kRed);
      knots->SetMarkerSize(1.5);
      knots->Draw("s:u:v", "type==1", "same"); // knots

    }//End draw
    if (!ask()) {
        break;
      }
    
  }//end if seed */
  std::cout << "testfunction!" << std::endl;
  return 0;
} //END END TESTFUNCTION 2D
#endif
  

template <typename DataT, bool isConsistentT>
int SplineBase<DataT, isConsistentT>::test(const bool draw, const bool drawDataPoints)
{
  // Test method

  using namespace std;

  constexpr int nDimX = 2;
  constexpr int nDimY = 2;
  constexpr int Fdegree = 4;

  double xMin[nDimX];
  double xMax[nDimX];
  int nKnots[nDimX];
  int* knotsU[nDimX];
  int nAxiliaryDatapoints[nDimX];

  for (int i = 0; i < nDimX; i++) {
    xMin[i] = 0.;
    xMax[i] = 1.;
    nKnots[i] = 4;
    knotsU[i] = new int[nKnots[i]];
    nAxiliaryDatapoints[i] = 4;
  }

  // Function F
  const int nTerms1D = 2 * (Fdegree + 1);
  int nFcoeff = nDimY;
  for (int i = 0; i < nDimX; i++) {
    nFcoeff *= nTerms1D;
  }

  double Fcoeff[nFcoeff];

  auto F = [&](const double x[nDimX], double f[nDimY]) {
    double a[nFcoeff];
    a[0] = 1;
    int na = 1;
    for (int d = 0; d < nDimX; d++) {
      double b[nFcoeff];
      int nb = 0;
      double t = (x[d] - xMin[d]) * TMath::Pi() / (xMax[d] - xMin[d]);
      for (int i = 0; i < nTerms1D; i++) {
        double c = (i % 2) ? cos((i / 2) * t) : cos((i / 2) * t);
        for (int j = 0; j < na; j++) {
          b[nb++] = c * a[j];
          assert(nb <= nFcoeff);
        }
      }
      na = nb;
      for (int i = 0; i < nb; i++) {
        a[i] = b[i];
      }
    }

    double* c = Fcoeff;
    for (int dim = 0; dim < nDimY; dim++) {
      f[dim] = 0;
      for (int i = 0; i < na; i++) {
        f[dim] += a[i] * (*c++);
      }
    }
  };

  auto F2D = [&](double x1, double x2, double f[nDimY]) {
    double x[2] = {x1, x2};
    F(x, f);
  };

  for (int seed = 1; seed < 10; seed++) {

    gRandom->SetSeed(seed);

    // getting the coefficents filled randomly
    for (int i = 0; i < nFcoeff; i++) {
      Fcoeff[i] = gRandom->Uniform(-1, 1);
    }

    for (int i = 0; i < nDimX; i++) {
      knotsU[i][0] = 0;
      for (int j = 1; j < nKnots[i]; j++) {
        knotsU[i][j] = j * 4; //+ int(gRandom->Integer(3)) - 1;
      }
    }

    Spline<double, nDimX, nDimY> spline(nKnots, knotsU);
    Spline2D<double, nDimY> spline2D(nKnots[0], knotsU[0], nKnots[1], knotsU[1]);

    cout << "mark 1" << std::endl;
    spline.approximateFunction(xMin, xMax, F, nAxiliaryDatapoints);
    spline2D.approximateFunction(xMin[0], xMax[0], xMin[1], xMax[1],
                                 F2D, nAxiliaryDatapoints[0], nAxiliaryDatapoints[0]);
    cout << "mark 2" << std::endl;

    long double statDf = 0;
    long double statDf2D = 0;

    long double statN = 0;

    double x[nDimX];
    for (int i = 0; i < nDimX; i++) {
      x[i] = xMin[i];
    }
    do {
      double xf[nDimX];
      double s[nDimY];
      double s2D[nDimY];
      double f[nDimY];
      for (int i = 0; i < nDimX; i++) {
        xf[i] = x[i];
      }
      F(x, f);
      spline.interpolate(xf, s);
      spline2D.interpolate(xf[0], xf[1], s2D);

      for (int dim = 0; dim < nDimY; dim++) {
        statDf += (s[dim] - f[dim]) * (s[dim] - f[dim]);
        statDf2D += (s2D[dim] - f[dim]) * (s2D[dim] - f[dim]);
        statN++;
      }
      int dim = 0;
      for (; dim < nDimX; dim++) {
        x[dim] += 0.01;
        if (x[dim] <= xMax[dim]) {
          break;
        }
        x[dim] = xMin[dim];
      }
      if (dim >= nDimX) {
        break;
      }
    } while (1);

    cout << "\n std dev for SplineND   : " << sqrt(statDf / statN) << std::endl;
    cout << "\n std dev for Spline2D   : " << sqrt(statDf2D / statN) << std::endl;

  } // seed

  for (int i = 0; i < nDimX; i++) {
    delete[] knotsU[i];
  }

  return 0;
}



#endif // GPUCA_GPUCODE

template class GPUCA_NAMESPACE::gpu::SplineBase<float, false>;
template class GPUCA_NAMESPACE::gpu::SplineBase<float, true>;
template class GPUCA_NAMESPACE::gpu::SplineBase<double, false>;
template class GPUCA_NAMESPACE::gpu::SplineBase<double, true>;

template class GPUCA_NAMESPACE::gpu::Spline<float, 2, 3, true>;
template class GPUCA_NAMESPACE::gpu::Spline<float, 2, 2, true>;
