/*

   // works only with ROOT >= 6

   alienv load ROOT/latest-root6
   alienv load Vc/latest

   root -l

   .x loadlibs.C
   .x IrregularSpline1DTest.C++
 */
int CompactSplineIrregular1DTest()
{
  return 0;
}
#ifdef XXX

#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "TFile.h"
#include "TRandom.h"
#include "TNtuple.h"
#include "Riostream.h"
#include "TSystem.h"
#include "TH1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TLine.h"
#include "GPU/CompactSplineIrregular1D.h"
#include "GPU/IrregularSpline1D.h"
#include "GPU/CompactSplineHelper.h"

const int funcN = 10;
static double funcC[2 * funcN + 2];

float F(float u)
{
  double f = 0; //funcC[0]/2;
  double uu = u * TMath::Pi();
  for (int i = 1; i <= funcN; i++) {
    f += funcC[2 * i] * TMath::Cos(i * uu) + funcC[2 * i + 1] * TMath::Sin(i * uu);
  }
  return f;
}

int CompactSplineIrregular1DTest()
{

  const int nAxiliaryPoints = 10;

  using namespace GPUCA_NAMESPACE::gpu;

  cout << "Test interpolation.." << endl;

  gRandom->SetSeed(0);

  for (int seed = 14;; seed++) {

    //seed = gRandom->Integer(100000); // 605

    gRandom->SetSeed(seed);
    cout << "Random seed: " << seed << " " << gRandom->GetSeed() << endl;

    for (int i = 0; i <= funcN; i++) {
      funcC[i] = gRandom->Uniform(-1, 1);
    }

    const int initNKnotsU = 5;
    int nAxisBinsU = 100;

    float knotsU[initNKnotsU];
    {
      knotsU[0] = 0;
      float du = 1. / (initNKnotsU - 1);

      for (int i = 1; i < initNKnotsU; i++) {
        knotsU[i] = i * du; // + gRandom->Uniform(-du / 3, du / 3);
      }
      knotsU[initNKnotsU - 1] = 1;
    }

    CompactSplineHelper helper;

    CompactSplineIrregular1D spline;
    spline.construct(initNKnotsU, knotsU, nAxisBinsU);
    std::unique_ptr<float[]> data = helper.create(spline, F, nAxiliaryPoints);

    int nKnotsTot = spline.getNumberOfKnots();
    cout << "Knots: initial " << initNKnotsU << ", created " << nKnotsTot << endl;
    for (int i = 0; i < nKnotsTot; i++) {
      cout << "knot " << i << ": " << spline.getKnot(i).u << endl;
    }

    CompactSplineIrregular1D splineClassic;
    splineClassic.construct(initNKnotsU, knotsU, nAxisBinsU);
    std::unique_ptr<float[]> dataClassic = helper.createClassical(splineClassic, F);

    IrregularSpline1D splineLocal;
    splineLocal.construct(initNKnotsU, knotsU, nAxisBinsU + 1);
    std::unique_ptr<float[]> dataLocal(new float[splineLocal.getNumberOfKnots()]);
    {
      for (int i = 0; i < splineLocal.getNumberOfKnots(); i++) {
        dataLocal[i] = F(splineLocal.getKnot(i).u);
      }
      splineLocal.correctEdges(dataLocal.get());
    }

    IrregularSpline1D splineLocal2N;
    std::unique_ptr<float[]> dataLocal2N(nullptr);
    {
      float knotsUtmp[nKnotsTot * 2];
      float du = 1. / (nKnotsTot * 2 - 1);
      for (int i = 0; i < nKnotsTot * 2; i++) {
        knotsUtmp[i] = i * du;
      }
      knotsUtmp[0] = 0.;
      knotsUtmp[2 * nKnotsTot - 1] = 1.;
      splineLocal2N.construct(nKnotsTot * 2, knotsUtmp, nAxisBinsU * 2 + 1);

      int n = splineLocal2N.getNumberOfKnots();
      dataLocal2N.reset(new float[n]); // corrected data

      for (int i = 0; i < n; i++) {
        dataLocal2N[i] = F(splineLocal2N.getKnot(i).u);
      }
      splineLocal2N.correctEdges(dataLocal2N.get());
    }

    spline.print();
    splineLocal.print();
    splineLocal2N.print();

    const CompactSplineIrregular1D& gridU = spline;
    int nu = gridU.getNumberOfKnots();

    TCanvas* canv = new TCanvas("cQA", "CompactSplineIrregular1D  QA", 2000, 1000);
    canv->Draw();

    TH1F* qaX = new TH1F("qaX", "qaX [um]", 1000, -1000., 1000.);

    TNtuple* knots = new TNtuple("knots", "knots", "type:u:f");

    double diff = 0;
    for (int i = 0; i < nu; i++) {
      double u = gridU.getKnot(i).u;
      double f0 = F(u);
      double fs = spline.getSpline((const float*)data.get(), u);
      diff += (fs - f0) * (fs - f0);
      knots->Fill(1, u, fs);

      if (i < nu - 1) {
        double u1 = gridU.getKnot(i + 1).u;
        int nax = nAxiliaryPoints;
        double du = (u1 - u) / (nax + 1);
        for (int j = 0; j < nax; j++) {
          double uu = u + du * (j + 1);
          double ff = spline.getSpline((const float*)data.get(), uu);
          knots->Fill(2, uu, ff);
        }
      }
    }

    cout << "mean diff at knots: " << sqrt(diff) / nu << endl;

    for (int i = 0; i < splineLocal2N.getNumberOfKnots(); i++) {
      double u = splineLocal2N.getKnot(i).u;
      double fs = splineLocal2N.getSpline((const float*)dataLocal2N.get(), u);
      knots->Fill(3, u, fs);
    }

    TNtuple* nt = new TNtuple("nt", "nt", "u:f0:fComp:fClass:fLocal:fLocal2N");

    float stepu = 1.e-4;
    int nSteps = (int)(1. / stepu + 1);

    double statDfLocal = 0;
    double statDfComp = 0;
    double statDfClass = 0;
    double statDfLocal2N = 0;
    int statN = 0;
    for (float u = 0; u < 1. + stepu; u += stepu) {
      double f0 = F(u);
      double fComp = spline.getSpline((const float*)data.get(), u);
      double fClass = splineClassic.getSpline((const float*)dataClassic.get(), u);
      double fLocal = splineLocal.getSpline((const float*)dataLocal.get(), u);
      double fLocal2N = splineLocal2N.getSpline((const float*)dataLocal2N.get(), u);

      statDfComp += (fComp - f0) * (fComp - f0);
      statDfClass += (fClass - f0) * (fClass - f0);
      statDfLocal += (fLocal - f0) * (fLocal - f0);
      statDfLocal2N += (fLocal2N - f0) * (fLocal2N - f0);
      statN++;
      qaX->Fill(1.e4 * (fComp - f0));
      nt->Fill(u, f0, fComp, fClass, fLocal, fLocal2N);
    }

    cout << "\n"
         << std::endl;
    cout << "std dev Classical : " << sqrt(statDfClass / statN) << std::endl;
    cout << "std dev Local     : " << sqrt(statDfLocal / statN) << std::endl;
    cout << "std dev Local 2N  : " << sqrt(statDfLocal2N / statN) << std::endl;
    cout << "std dev Compact   : " << sqrt(statDfComp / statN) << std::endl;

    /*
      canv->cd(1);
      qaX->Draw();
      canv->cd(2);
    */

    //nt->SetMarkerColor(kBlack);
    //nt->Draw("f0:u","","");

    nt->SetMarkerColor(kGray);
    nt->SetMarkerStyle(8);
    nt->SetMarkerSize(2.);
    nt->Draw("f0:u", "", "P");

    TH1* htemp = (TH1*)gPad->GetPrimitive("htemp");
    htemp->SetTitle("Splines of the same size");

    knots->SetMarkerStyle(8);
    knots->SetMarkerSize(1.5);

    nt->SetMarkerSize(.5);
    nt->SetMarkerColor(kBlack);
    nt->Draw("fClass:u", "", "P,same");

    /*
      nt->SetMarkerColor(kCyan);
      nt->Draw("fLocal:u", "", "P,same");
    */

    nt->SetMarkerColor(kBlue);
    nt->Draw("fLocal2N:u", "", "P,same");

    nt->SetMarkerColor(kRed);
    nt->Draw("fComp:u", "", "P,same");

    knots->SetMarkerSize(2.);
    knots->SetMarkerColor(kBlue);
    knots->Draw("f:u", "type==3", "same");

    knots->SetMarkerColor(kBlack);
    knots->SetMarkerSize(2.5);
    knots->Draw("f:u", "type==1", "same");

    knots->SetMarkerColor(kRed);
    knots->SetMarkerSize(1.5);
    knots->Draw("f:u", "type==1", "same");

    knots->SetMarkerColor(kBlack);
    knots->SetMarkerSize(1.);
    knots->Draw("f:u", "type==2", "same");

    auto legend = new TLegend(0.1, 0.85, 0.3, 0.95);
    //legend->SetHeader("Splines of the same size:","C"); // option "C" allows to center the header
    TLine l1;
    l1.SetLineWidth(4);
    l1.SetLineColor(kBlack);
    legend->AddEntry(&l1, "classical (N knots + N slopes)", "L");
    TLine l2(l1);
    l2.SetLineColor(kBlue);
    legend->AddEntry(&l2, "local (2N knots)", "L");
    TLine l3(l1);
    l3.SetLineColor(kRed);
    legend->AddEntry(&l3, "compact (N knots + N slopes)", "L");
    legend->Draw();

    canv->Update();

    cout << "\nRandom seed: " << seed << " " << gRandom->GetSeed() << endl;
    cout << "type 'q' to exit" << endl;
    if (getchar() == 'q')
      break;
  }

  return 0;
}

#endif
#endif