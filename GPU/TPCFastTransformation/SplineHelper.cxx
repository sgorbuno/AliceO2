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
/// \brief Implementation of SplineHelper class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE)

#include "SplineHelper.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompBK.h"

#include <vector>
#include "TRandom.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "TFile.h"
#include "GPUCommonMath.h"
#include <iostream>

using namespace GPUCA_NAMESPACE::gpu;

template <typename DataT>
SplineHelper<DataT>::SplineHelper() : mError(), mXdimensions(0), mFdimensions(0), mNumberOfDataPoints(0), mHelpers()
{
}

template <typename DataT>
int SplineHelper<DataT>::storeError(int code, const char* msg)
{
  mError = msg;
  return code;
}

////////////////
// pointstoarray
// HILFSFUNKTION,
template <typename DataT>
int SplineHelper<DataT>::pointstoarray(const int indices[], const int numbers[], int dim)
{

  int result = 0;
  int factor = 1;
  for (int i = 0; i < dim; i++) {
    result += indices[i] * factor;
    factor *= numbers[i];
  }
  return result;
}

////////////////
//arraytopoints
// HILFSFUNKTION
template <typename DataT>
int SplineHelper<DataT>::arraytopoints(int point, int result[], const int numbers[], int dim)
{

  if (point == 0) {
    for (int i = 0; i < dim; i++) {
      result[i] = 0;
    }
  } else {
    int divisor = 1;
    int modoperand = 1;
    for (int i = 0; i < dim; i++) {
      modoperand *= numbers[i];
      result[i] = (int)((point % modoperand) / divisor);
      divisor *= numbers[i];
    }
  }
  return 0;
}

template <typename DataT>
void SplineHelper<DataT>::approximateFunction(
  DataT* Fparameters, const double xMin[/* mXdimensions */], const double xMax[/* mXdimensions */],
  std::function<void(const double x[/* mXdimensions */], double f[/* mFdimensions */])> F) const
{
  /// Create best-fit spline parameters for a given input function F
  /// output in Fparameter
  // TODO: implement
  // MY VERSION
  //std::cout << "approximateFunction(Fparameters, xMin[],xMax[],F) :" << std::endl;
  double scaleX[mXdimensions];
  for (int i = 0; i < mXdimensions; i++) {
    scaleX[i] = (xMax[i] - xMin[i]) / ((double)(mHelpers[i].getSpline().getUmax()));
  }

  // calculate F-Values at all datapoints:
  int nrOfAllPoints = getNumberOfDataPoints();
  std::vector<double> dataPointF(nrOfAllPoints * mFdimensions);

  int nrOfPoints[mXdimensions];
  for (int i = 0; i < mXdimensions; i++) {
    nrOfPoints[i] = mHelpers[i].getNumberOfDataPoints();
  }
  double x[mXdimensions];
  for (int d = 0; d < nrOfAllPoints; d++) { // for all DataPoints

    int indices[mXdimensions];
    int modoperand = 1;
    int divisor = 1;

    // get the DataPoint index
    for (int i = 0; i < mXdimensions; i++) {
      modoperand *= nrOfPoints[i];

      indices[i] = (int)((d % modoperand) / divisor);
      divisor *= nrOfPoints[i];
      // get the respecting u-values:
      x[i] = xMin[i] + mHelpers[i].getDataPoint(indices[i]).u * scaleX[i];
    }

    for (int j = 0; j < mXdimensions; j++) {
      F(x, &dataPointF[d * mFdimensions]);
    }

  } // end for all DataPoints d
  //END MY VERSION

  //std::vector<DataT> dataPointF(getNumberOfDataPoints() * mFdimensions);
  //DUMYY VERSION Commented out
  /* for (int i = 0; i < getNumberOfDataPoints() * mFdimensions; i++) {
    dataPointF[i] = 1.; 
  } */
  /*
  double scaleX1 = (x1Max - x1Min) / ((double)mHelperU1.getSpline().getUmax());
  double scaleX2 = (x2Max - x2Min) / ((double)mHelperU2.getSpline().getUmax());

  for (int iv = 0; iv < getNumberOfDataPointsU2(); iv++) {
    DataT x2 = x2Min + mHelperU2.getDataPoint(iv).u * scaleX2;
    for (int iu = 0; iu < getNumberOfDataPointsU1(); iu++) {
      DataT x1 = x1Min + mHelperU1.getDataPoint(iu).u * scaleX1;
      F(x1, x2, &dataPointF[(iv * getNumberOfDataPointsU1() + iu) * mFdimensions]);
    }
  }
  */
  approximateFunction(Fparameters, dataPointF.data());
}

template <typename DataT>
void SplineHelper<DataT>::approximateFunctionBatch(
  DataT* Fparameters, const double xMin[], const double xMax[],
  std::function<void(const std::vector<double> x[], std::vector<double> f[/*mFdimensions*/])> F,
  unsigned int batchsize) const
{
  /// Create best-fit spline parameters for a given input function F.
  /// F calculates values for a batch of points.
  /// output in Fparameters

  // TODO: implement later

  std::vector<double> dataPointF(getNumberOfDataPoints() * mFdimensions);
  for (int i = 0; i < getNumberOfDataPoints() * mFdimensions; i++) {
    dataPointF[i] = 0.;
  }

  approximateFunction(Fparameters, dataPointF.data());
}

template <typename DataT>
void SplineHelper<DataT>::approximateFunction(
  DataT* Fparameters, const double DataPointF[/*getNumberOfDataPoints() x nFdim*/]) const
{
  /// approximate a function given as an array of values at data points

  int numberOfKnots[mXdimensions]; //getting number of Knots for all dimensions into one array
  for (int i = 0; i < mXdimensions; i++) {
    numberOfKnots[i] = mHelpers[i].getSpline().getNumberOfKnots();
  }

  int numberOfDataPoints[mXdimensions]; // getting number of datapoints (incl knots) in all dimensions into one array
  for (int i = 0; i < mXdimensions; i++) {
    numberOfDataPoints[i] = mHelpers[i].getNumberOfDataPoints();
  }

  int numberOfAllKnots = 1; // getting Number of all knots for the entire spline
  for (int i = 0; i < mXdimensions; i++) {
    numberOfAllKnots *= numberOfKnots[i];
  }
  // TO BE REMOVED (TEST-OUTPUT):
  std::cout << "total number of knots: " << numberOfAllKnots << ", " << std::endl;

  int numberOfAllDataPoints = 1; // getting Number of all Datapoints for the entire spline
  for (int i = 0; i < mXdimensions; i++) {
    numberOfAllDataPoints *= numberOfDataPoints[i];
    //std::cout << mHelpers[0].getNumberOfDataPoints()<<std::endl;
  }

  // TO BE REMOVED TEST:
  //std::cout << "total number of DataPoints (including knots): " <<  numberOfAllDataPoints << ", "<< std::endl;

  int numberOfParameterTypes = (int)(pow(2.0, mXdimensions)); //number of Parameters per Knot

  // TO BE REMOVED TEST:
  //std::cout << "number of paramtertypes per knot : " <<  numberOfParameterTypes << ", "<< std::endl;

  std::unique_ptr<double[]> allParameters[numberOfParameterTypes]; //Array for the different parametertypes s, s'u, s'v, s''uv,...
  for (int i = 0; i < numberOfParameterTypes; i++) {
    allParameters[i] = std::unique_ptr<double[]>(new double[numberOfAllDataPoints * mFdimensions]); //To-Do:Fdim!!
  }
  //filling allParameters[0] and FParameters with s:
  for (int i = 0; i < numberOfAllDataPoints; i++) {
    for (int f = 0; f < mFdimensions; f++) {                                     // for all f-dimensions
      allParameters[0][i * mFdimensions + f] = DataPointF[i * mFdimensions + f]; // TO DO - Just get the pointer adress there PLEASE!
    }
    int p0indices[mXdimensions];
    arraytopoints(i, p0indices, numberOfDataPoints, mXdimensions);
    bool isKnot = 1;
    for (int j = 0; j < mXdimensions; j++) { // is the current datapoint a knot?
      if (!mHelpers[j].getDataPoint(p0indices[j]).isKnot) {
        isKnot = 0;
        break;
      }
    }
    if (isKnot) {
      int knotindices[mXdimensions];
      for (int j = 0; j < mXdimensions; j++) { // calculate KNotindices for all dimensions
        //WORKAROUND Getting Knotindices:
        knotindices[j] = p0indices[j] / ((numberOfDataPoints[j] - 1) / (numberOfKnots[j] - 1));
        //knotindices[j] = mHelpers[j].getDataPoint(p0indices[j]).iKnot; //in der Annahme der wert ist ein Knotenindex und falls der datapoint ein knoten ist, gibt er seinen eigenen knotenindex zurück
      }
      // get the knotindexvalue for FParameters:
      int knotind = pointstoarray(knotindices, numberOfKnots, mXdimensions);

      for (int f = 0; f < mFdimensions; f++) {                                                               // for all f-dimensions get function values into Fparameters
        Fparameters[knotind * numberOfParameterTypes * mFdimensions + f] = DataPointF[i * mFdimensions + f]; ///write derivatives in FParameters
      }
    } // end if isKnot
  }   //end i (filling DataPointF Values into allParameters[0] and FParameters)
  //now: allParameters[0] = dataPointF;

  //Array for input DataPointF-values for Spline1D::approximateFunctionGradually(...);
  std::unique_ptr<double[]> dataPointF1D[mXdimensions];
  for (int i = 0; i < mXdimensions; i++) {
    dataPointF1D[i] = std::unique_ptr<double[]>(new double[numberOfDataPoints[i] * mFdimensions]); // To-Do:Fdim!! For s and derivetives at all knots.
  }
  //Array to be filled by Spline1D::approximateFunctionGradually(...);
  std::unique_ptr<DataT[]> par[mXdimensions];
  std::unique_ptr<double[]> parD[mXdimensions];

  for (int i = 0; i < mXdimensions; i++) {
    par[i] = std::unique_ptr<DataT[]>(new DataT[numberOfKnots[i] * mFdimensions * 2]);
    parD[i] = std::unique_ptr<double[]>(new double[numberOfKnots[i] * mFdimensions * 2]);
  }

  //std::cout << "NumberOfParameters: " <<  mNumberOfParameters <<std::endl;

  //STARTING MAIN-LOOP, for all Parametertypes:
  for (int p = 1; p < numberOfParameterTypes; p++) { // p = 1!! Wir kriegen s (p0) durch approximateFunction()oben
    int dimension;                                   // find the dimension for approximation
    for (int i = (int)(log2f((float)p)); i >= 0; i--) {
      if (p % (int)(pow(2.0, i)) == 0) {
        dimension = i;
        break;
      }
    }

    int currentDataPointF = p - (int)(pow(2.0, dimension));
    //std::cout << std::endl << "p:" << p << ", dim of approximation: " << dimension << ", based on: " << currentDataPointF << std::endl;

    int nrOf1DSplines = (numberOfAllDataPoints / numberOfDataPoints[dimension]); // number of Splines for Parametertyp p in direction dim
    //std::cout << "nr of splines: " << nrOf1DSplines;

    // getting the numbers of Datapoints for all dimension eccept the dimension of interpolation
    int currentNumbers[mXdimensions - 1];
    for (int i = 0; i < dimension; i++) {
      currentNumbers[i] = numberOfDataPoints[i];
    }
    for (int i = dimension; i < mXdimensions - 1; i++) {
      currentNumbers[i] = numberOfDataPoints[i + 1];
    }
    std::cout << " current numbers: ";
    for (int i = 0; i < mXdimensions - 1; i++) {
      //std::cout << currentNumbers[i] << ",";
    }
    //std::cout << std::endl;

    //// for all Splines in current dimension:
    for (int s = 0; s < nrOf1DSplines; s++) {
      int indices[mXdimensions - 1];
      arraytopoints(s, indices, currentNumbers, mXdimensions - 1);
      int startpoint[mXdimensions]; //startpoint for the current 1DSpline
      for (int i = 0; i < dimension; i++) {
        startpoint[i] = indices[i];
      }
      startpoint[dimension] = 0;
      for (int i = dimension + 1; i < mXdimensions; i++) {
        startpoint[i] = indices[i - 1];
      }
      // NOW WE HAVE THE DATAPOINTINDICES OF THE CURRENT STARTPOINT IN startpoint-Array.
      int startdatapoint = pointstoarray(startpoint, numberOfDataPoints, mXdimensions);
      int distance = 1; // distance to the next dataPoint in the array for the current dimension
      for (int i = 0; i < dimension; i++) {
        distance *= numberOfDataPoints[i];
      }
      distance *= mFdimensions;

      for (int i = 0; i < numberOfDataPoints[dimension]; i++) { //Fill the dataPointF1D-Array
        for (int f = 0; f < mFdimensions; f++) {
          dataPointF1D[dimension][i * mFdimensions + f] = allParameters[currentDataPointF][startdatapoint * mFdimensions + (i * distance + f)]; // uiuiui index kuddelmuddel???!!
        }
      }
      mHelpers[dimension].approximateFunction(par[dimension].get(), dataPointF1D[dimension].get());
      for (int i = 0; i < numberOfKnots[dimension] * mFdimensions * 2; i++) {
        parD[dimension][i] = par[dimension][i];
      }
      // now we have all s and s' values in par[dimension]

      int redistributionindex[mXdimensions];
      for (int i = 0; i < mXdimensions; i++) {
        redistributionindex[i] = startpoint[i];
      }
      //redistributing the derivatives at dimension-Knots into array p
      for (int i = 0; i < numberOfKnots[dimension]; i++) {                        //for all dimension-Knots
        redistributionindex[dimension] = mHelpers[dimension].getKnotDataPoint(i); //find the indices
        int finalposition = pointstoarray(redistributionindex, numberOfDataPoints, mXdimensions);

        for (int f = 0; f < mFdimensions; f++) {
          allParameters[p][finalposition * mFdimensions + f] = par[dimension][2 * i * mFdimensions + mFdimensions + f];
        }

        bool isKnot = 1;
        for (int j = 0; j < mXdimensions; j++) { //is dataPoint a knot?
          if (!mHelpers[j].getDataPoint(redistributionindex[j]).isKnot) {
            isKnot = 0;
            break;
          } //noch mal checken!! Das muss noch anders!!
        }

        if (isKnot) { // for all knots
          int knotindices[mXdimensions];

          for (int j = 0; j < mXdimensions; j++) { // calculate Knotindices for all dimensions
            knotindices[j] = redistributionindex[j] / ((numberOfDataPoints[j] - 1) / (numberOfKnots[j] - 1));
            //knotindices[j] = mHelpers[j].getDataPoint(redistributionindex[j]).iKnot; //in der Annahme der wert ist ein Knotenindex und falls der datapoint ein knoten ist, gibt er seinen eigenen knotenindex zurück
          }
          // get the knotindexvalue for FParameters:
          int knotind = pointstoarray(knotindices, numberOfKnots, mXdimensions);
          for (int f = 0; f < mFdimensions; f++)
            Fparameters[knotind * numberOfParameterTypes * mFdimensions + p * mFdimensions + f] = par[dimension][2 * i * mFdimensions + mFdimensions + f]; ///write derivatives in FParameters
        }
      } // end for all fknots (for redistribution)

      // recalculation:
      for (int i = 0; i < numberOfDataPoints[dimension]; i++) { // this is somehow still redundant// TO DO: ONLY PART OF approximateFunction WHERE NDIM is considerd!!
        redistributionindex[dimension] = i;                     // getting current datapointindices
        bool isKnot = 1;                                        // check is current datapoint a knot?
        for (int j = 0; j < mXdimensions; j++) {
          if (!mHelpers[j].getDataPoint(redistributionindex[j]).isKnot) {
            isKnot = 0;
            break;
          }
        }
        double splineF[mFdimensions];
        double u = mHelpers[dimension].getDataPoint(i).u;
        mHelpers[dimension].getSpline().interpolateU(mFdimensions, parD[dimension].get(), u, splineF); //recalculate at all datapoints of dimension
        for (int dim = 0; dim < mFdimensions; dim++) {                                                //writing it in allParameters
          //std::cout<<allParameters [p-(int)(pow(2.0, dimension))] [(int)(startdatapoint*mFdimensions + i*distance + dim)]<<", ";
          allParameters[p - (int)(pow(2.0, dimension))][(int)(startdatapoint * mFdimensions + i * distance + dim)] = splineF[dim]; //write it in the array.
          //std::cout<<allParameters [p-(int)(pow(2.0, dimension))] [(int)(startdatapoint*mFdimensions + i*distance + dim)]<<",   ";
        }

        if (isKnot) {
          int knotindices[mXdimensions];

          for (int j = 0; j < mXdimensions; j++) { // calculate KNotindices for all dimensions
            knotindices[j] = redistributionindex[j] / ((numberOfDataPoints[j] - 1) / (numberOfKnots[j] - 1));
            //knotindices[j] = mHelpers[j].getDataPoint(redistributionindex[j]).iKnot; //in der Annahme der wert ist ein Knotenindex und falls der datapoint ein knoten ist, gibt er seinen eigenen knotenindex zurück
          }
          int currentknotarrayindex = pointstoarray(knotindices, numberOfKnots, mXdimensions);
          // getting the recalculated value into FParameters:
          for (int f = 0; f < mFdimensions; f++) {
            Fparameters[currentknotarrayindex * numberOfParameterTypes * mFdimensions + (p - (int)(pow(2.0, dimension))) * mFdimensions + f] = splineF[f];
          }
        } // end if isKnot
      }   // end recalculation
    }     //end of all1DSplines
  }       //end of for parametertypes
} //end of approxymateFunction MYVERSION!

//****
//* TESTFUNCTION 2D
template <typename DataT>
int SplineHelper<DataT>::test(const bool draw, const bool drawDataPoints)
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

template class GPUCA_NAMESPACE::gpu::SplineHelper<float>;
template class GPUCA_NAMESPACE::gpu::SplineHelper<double>;

#endif
