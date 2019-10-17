/*
   root -l -q CompactSplineDemo.C+
 */

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
#include "GPU/CompactSplineHelper1D.h"

const int funcN = 10;
static double funcC[2 * funcN + 2];

int nKnots = 4;

float F(float u)
{
  double uu = u * TMath::Pi() / (nKnots - 1);
  double f = 0; //funcC[0]/2;
  for (int i = 1; i <= funcN; i++) {
    f += funcC[2 * i] * TMath::Cos(i * uu) + funcC[2 * i + 1] * TMath::Sin(i * uu);
  }
  return f;
}

float F01(float u)
{
  return F(u * (nKnots - 1));
}

TCanvas* canv = new TCanvas("cQA", "Compact Spline Demo", 2000, 800);

bool doAskSteps = 1;

bool ask()
{
  canv->Update();
  if (doAskSteps)
    cout << "type 'q ' to exit, 's' to skip individual steps" << endl;
  else
    cout << "type 'q ' to exit, 's' to stop at individual steps" << endl;
  std::string str;
  std::getline(std::cin, str);
  if (str == "s")
    doAskSteps = !doAskSteps;
  return (str != "q" && str != ".q");
}

bool askStep()
{
  if (!doAskSteps)
    return 1;
  return ask();
}

int CompactSplineDemo()
{

  const int nAxiliaryPoints = 10;

  using namespace GPUCA_NAMESPACE::gpu;

  cout << "Test interpolation.." << endl;

  //TCanvas* canv = new TCanvas("cQA", "CompactSplineIrregular1D  QA", 2000, 1000);

  gRandom->SetSeed(0);

  for (int seed = 19;; seed++) {

    //seed = gRandom->Integer(100000); // 605

    gRandom->SetSeed(seed);
    cout << "Random seed: " << seed << " " << gRandom->GetSeed() << endl;

    for (int i = 0; i <= funcN; i++) {
      funcC[i] = gRandom->Uniform(-1, 1);
    }

    CompactSplineHelper1D helper;

    CompactSplineIrregular1D spline;
    spline.constructRegular(nKnots);
    std::unique_ptr<float[]> data = helper.constructSpline(spline, F, 0., spline.getUmax(), nAxiliaryPoints);

    CompactSplineIrregular1D splineClassic;
    splineClassic.constructRegular(nKnots);
    std::unique_ptr<float[]> dataClassic = helper.constructSplineClassical(splineClassic, F, 0., splineClassic.getUmax());

    IrregularSpline1D splineLocal;
    int nKnotsLocal = 2 * nKnots - 1;
    splineLocal.constructRegular(nKnotsLocal);

    std::unique_ptr<float[]> dataLocal(new float[nKnotsLocal]);
    for (int i = 0; i < nKnotsLocal; i++) {
      dataLocal[i] = F01(splineLocal.getKnot(i).u);
    }
    splineLocal.correctEdges(dataLocal.get());

    spline.print();
    splineLocal.print();

    canv->Draw();

    TH1F* qaX = new TH1F("qaX", "qaX [um]", 1000, -1000., 1000.);

    TNtuple* knots = new TNtuple("knots", "knots", "type:u:f");

    for (int i = 0; i < nKnots; i++) {
      double u = splineClassic.getKnot(i).u;
      double fs = splineClassic.getSpline((const float*)dataClassic.get(), u);
      knots->Fill(1, u, fs);
    }

    for (int i = 0; i < nKnots; i++) {
      double u = spline.getKnot(i).u;
      double fs = spline.getSpline((const float*)data.get(), u);
      knots->Fill(2, u, fs);
      if (i < nKnots - 1) {
        double u1 = spline.getKnot(i + 1).u;
        int nax = nAxiliaryPoints;
        double du = (u1 - u) / (nax + 1);
        for (int j = 0; j < nax; j++) {
          double uu = u + du * (j + 1);
          double ff = spline.getSpline((const float*)data.get(), uu);
          knots->Fill(3, uu, ff);
        }
      }
    }

    for (int i = 0; i < splineLocal.getNumberOfKnots(); i++) {
      double u = splineLocal.getKnot(i).u;
      double fs = splineLocal.getSpline((const float*)dataLocal.get(), u);
      knots->Fill(4, u * (nKnots - 1), fs);
    }

    TNtuple* nt = new TNtuple("nt", "nt", "u:f0:fComp:fClass:fLocal");

    float stepS = 1.e-4;
    int nSteps = (int)(1. / stepS + 1);

    double statDfComp = 0;
    double statDfClass = 0;
    double statDfLocal = 0;

    int statN = 0;
    for (float s = 0; s < 1. + stepS; s += stepS) {
      double u = s * (nKnots - 1);
      double f0 = F(u);
      double fComp = spline.getSpline((const float*)data.get(), u);
      double fClass = splineClassic.getSpline((const float*)dataClassic.get(), u);
      double fLocal = splineLocal.getSpline((const float*)dataLocal.get(), s);

      statDfComp += (fComp - f0) * (fComp - f0);
      statDfClass += (fClass - f0) * (fClass - f0);
      statDfLocal += (fLocal - f0) * (fLocal - f0);
      statN++;
      qaX->Fill(1.e4 * (fComp - f0));
      nt->Fill(u, f0, fComp, fClass, fLocal);
    }

    cout << "\n"
         << std::endl;
    cout << "\nRandom seed: " << seed << " " << gRandom->GetSeed() << endl;
    cout << "std dev Classical : " << sqrt(statDfClass / statN) << std::endl;
    cout << "std dev Local     : " << sqrt(statDfLocal / statN) << std::endl;
    cout << "std dev Compact   : " << sqrt(statDfComp / statN) << std::endl;

    /*
      canv->cd(1);
      qaX->Draw();
      canv->cd(2);
    */

    //nt->SetMarkerColor(kBlack);
    //nt->Draw("f0:u","","");

    auto legend = new TLegend(0.1, 0.82, 0.3, 0.95);
    //legend->SetHeader("Splines of the same size:","C"); // option "C" allows to center the header

    nt->SetMarkerColor(kGray);
    nt->SetMarkerStyle(8);
    nt->SetMarkerSize(2.);
    nt->Draw("f0:u", "", "P");

    TH1* htemp = (TH1*)gPad->GetPrimitive("htemp");
    htemp->SetTitle("Splines of the same size");

    TLine* l0 = new TLine();
    l0->SetLineWidth(10);
    l0->SetLineColor(kGray);
    legend->AddEntry(l0, "input function", "L");
    legend->Draw();

    knots->SetMarkerStyle(8);
    knots->SetMarkerSize(1.5);

    if (!askStep())
      break;

    nt->SetMarkerSize(.5);
    nt->SetMarkerColor(kGreen + 2);
    nt->Draw("fClass:u", "", "P,same");

    knots->SetMarkerColor(kGreen + 2);
    knots->SetMarkerSize(3.5);
    knots->Draw("f:u", "type==1", "same"); // classical
    TLine* l1 = new TLine();
    l1->SetLineWidth(4);
    l1->SetLineColor(kGreen + 2);
    legend->AddEntry(l1, "classical (N knots + N slopes)", "L");
    legend->Draw();

    if (!askStep())
      break;

    nt->SetMarkerColor(kBlue);
    nt->Draw("fLocal:u", "", "P,same");

    knots->SetMarkerSize(2.5);
    knots->SetMarkerColor(kBlue);
    knots->Draw("f:u", "type==4", "same"); // local
    TLine* l2 = new TLine(*l1);
    l2->SetLineColor(kBlue);
    legend->AddEntry(l2, "local (2N knots)", "L");
    legend->Draw();
    if (!askStep())
      break;

    nt->SetMarkerColor(kRed);
    nt->Draw("fComp:u", "", "P,same");

    knots->SetMarkerColor(kRed);
    knots->SetMarkerSize(1.5);
    knots->Draw("f:u", "type==2", "same"); // compact
    TLine* l3 = new TLine(*l1);
    l3->SetLineColor(kRed);
    legend->AddEntry(l3, "compact (N knots + N slopes)", "L");
    legend->Draw();

    if (!askStep())
      break;

    knots->SetMarkerColor(kBlack);
    knots->SetMarkerSize(1.);
    knots->Draw("f:u", "type==3", "same"); // compact, axiliary points
    TMarker* l4 = new TMarker;
    l4->SetMarkerStyle(8);
    l4->SetMarkerSize(1.);
    l4->SetMarkerColor(kBlack);
    legend->AddEntry(l4, "construction points", "P");
    legend->Draw();

    if (!ask())
      break;
    delete legend;
  }

  return 0;
}

#endif