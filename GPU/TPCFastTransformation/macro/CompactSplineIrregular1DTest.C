/*
  root -l -q IrregularSpline1DTest.C+
 */

#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "TRandom.h"
#include "Riostream.h"
#include "TMath.h"
#include "GPU/CompactSplineIrregular1D.h"
#include "GPU/CompactSplineHelper.h"
#include "TCanvas.h"
#include "TNtuple.h"

const bool kDraw = 0;

const int funcN = 10;
static double funcC[2 * funcN + 2];

int nKnots = 4;
int uMax = nKnots * 3;

float F(float u)
{
  double uu = u * TMath::Pi() / uMax;
  double f = 0; //funcC[0]/2;
  for (int i = 1; i <= funcN; i++) {
    f += funcC[2 * i] * TMath::Cos(i * uu) + funcC[2 * i + 1] * TMath::Sin(i * uu);
  }
  return f;
}

TCanvas* canv = nullptr;

bool ask()
{
  if (!canv)
    return 0;
  canv->Update();
  cout << "type 'q ' to exit" << endl;
  std::string str;
  std::getline(std::cin, str);
  return (str != "q" && str != ".q");
}

int CompactSplineIrregular1DTest()
{
  const int nAxiliaryPoints = 5;

  using namespace GPUCA_NAMESPACE::gpu;

  cout << "Test 1D interpolation with the compact spline" << endl;

  int nTries = 100;

  if (kDraw) {
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
    cout << "std dev Compact   : " << sqrt(statDf / statN) << std::endl;

    if (kDraw) {
      TNtuple* nt = new TNtuple("nt", "nt", "u:f0:fSpline");
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

      TNtuple* knots = new TNtuple("knots", "knots", "type:u:f");
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
#else
int CompactSplineIrregular1DTest()
{
  std::cout << "\nPlease run this macro with precompilation: \"root -l -q IrregularSpline1DTest.C+\"\n"
            << std::endl;
  return 0;
}
#endif