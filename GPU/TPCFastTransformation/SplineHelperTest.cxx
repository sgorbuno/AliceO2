// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  Spline182.cxx
/// \brief Implementation of SplineHelperTest class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE)

#include "SplineHelperTest.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompBK.h"

using namespace GPUCA_NAMESPACE::gpu;

template <typename DataT>
SplineHelperTest<DataT>::SplineHelperTest() : mError(), mXdimensions(0), mFdimensions(0), mNumberOfDataPoints(0), mHelpers()
{
}

template <typename DataT>
int SplineHelperTest<DataT>::storeError(int code, const char* msg)
{
  mError = msg;
  return code;
}

////////////////
// pointstoarray
// HILFSFUNKTION,
template <typename DataT>
int SplineHelperTest<DataT>::pointstoarray ( const int indices[], const int numbers[], int dim){

  int result=0;
  int factor = 1;
  for (int i = 0; i<dim; i++){
    result += indices[i]*factor;
    factor *= numbers[i];
  }
  return result;
}

////////////////
//arraytopoints
// HILFSFUNKTION
template <typename DataT>
int SplineHelperTest<DataT>::arraytopoints ( int point, int result[], const int numbers[],  int dim){
  
  if (point == 0) {
      for (int i = 0; i < dim; i++){result[i]=0;}  
  }
  else{
    int divisor = 1;
    int modoperand = 1;
    for (int i = 0; i < dim; i++){
      modoperand *= numbers[i];
      result[i] = (int) ((point % modoperand)/divisor);
      divisor *= numbers[i];
    }
  }
  return 0;
}


template <typename DataT>
void SplineHelperTest<DataT>::approximateFunction(
  DataT* Fparameters, const DataT xMin[/* mXdimensions */], const DataT xMax[/* mXdimensions */],
  std::function<void(const DataT x[/* mXdimensions */], DataT f[/* mFdimensions */])> F) const
{
  /// Create best-fit spline parameters for a given input function F
  /// output in Fparameters

  // TODO: implement
  // MY VERSION
  std::cout << "approximateFunction(Fparameters, xMin[],xMax[],F) :" << std::endl;
  double scaleX[mXdimensions];
    for (int i = 0; i < mXdimensions; i++){
        scaleX[i] = (xMax[i] - xMin[i])/((double)(mHelpers[i].getSpline().getUmax()));
    }
    
    // calculate F-Values at all datapoints:
    int nrOfAllPoints = getNumberOfDataPoints();
    std::vector<DataT> dataPointF( nrOfAllPoints * mFdimensions);
    
    int nrOfPoints[mXdimensions];
    for (int i = 0; i < mXdimensions; i++){
      nrOfPoints[i] = mHelpers[i].getNumberOfDataPoints();
    }
    DataT x[mXdimensions];
    for (int d = 0; d < nrOfAllPoints; d++){    // for all DataPoints   

        int indices[mXdimensions];
        int modoperand = 1;
        int divisor = 1;
         
        // get the DataPoint index
        for (int i = 0; i < mXdimensions; i++){
            modoperand *= nrOfPoints[i];

            indices[i]= (int)((d % modoperand) /divisor);
            divisor *= nrOfPoints[i];
            // get the respecting u-values:      
            x[i]= xMin[i] + mHelpers[i].getDataPoint(indices[i]).u * scaleX[i];
        }
        
        for (int j = 0; j < mXdimensions; j++){
          F(x, &dataPointF[d * mFdimensions]);

        }
       
    }// end for all DataPoints d 
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
void SplineHelperTest<DataT>::approximateFunctionBatch(
  DataT* Fparameters, const DataT xMin[], const DataT xMax[],
  std::function<void(const std::vector<DataT> x[], std::vector<DataT> f[/*mFdimensions*/])> F,
  unsigned int batchsize) const
{
  /// Create best-fit spline parameters for a given input function F.
  /// F calculates values for a batch of points.
  /// output in Fparameters

  // TODO: implement later

  std::vector<DataT> dataPointF(getNumberOfDataPoints() * mFdimensions);
  for (int i = 0; i < getNumberOfDataPoints() * mFdimensions; i++) {
    dataPointF[i] = 0.;
  }

  approximateFunction(Fparameters, dataPointF.data());
}

template <typename DataT>
void SplineHelperTest<DataT>::approximateFunction(
  DataT* Fparameters, const DataT DataPointF[/*getNumberOfDataPoints() x nFdim*/]) const
{
  /// approximate a function given as an array of values at data points

  assert(mFdimensions == 1 );
  assert(mXdimensions == 2 );

  int numberOfKnots [mXdimensions]; //getting number of Knots for all dimensions into one array
  for (int i = 0; i < mXdimensions; i++){
      numberOfKnots[i]=mHelpers[i].getSpline().getNumberOfKnots();
  }
  
  int numberOfDataPoints[mXdimensions];// getting number of datapoints (incl knots) in all dimensions into one array
  for (int i = 0; i < mXdimensions; i++){
    numberOfDataPoints[i]=mHelpers[i].getNumberOfDataPoints();
  }

  int numberOfAllKnots = 1; // getting Number of all knots for the entire spline
  for (int i = 0; i < mXdimensions; i++){
    numberOfAllKnots *= numberOfKnots[i];
  }
  // TO BE REMOVED (TEST-OUTPUT):
  std::cout << "total number of knots: " << numberOfAllKnots << ", "<<std::endl;

  int numberOfAllDataPoints = 1; // getting Number of all Datapoints for the entire spline
  for (int i = 0; i < mXdimensions; i++){
    numberOfAllDataPoints *= numberOfDataPoints[i];
  }

  // TO BE REMOVED TEST:
  std::cout << "total number of DataPoints (including knots): " <<  numberOfAllDataPoints << ", "<< std::endl;
  
  int numberOfParameterTypes = (int)(pow(2.0, mXdimensions)); //number of Parameters per Knot
  
  // TO BE REMOVED TEST:
  std::cout << "number of paramtertypes per knot : " <<  numberOfParameterTypes << ", "<< std::endl;

 std::unique_ptr<DataT[]> allParameters[numberOfParameterTypes]; //Array for the different parametertypes s, s'u, s'v, s''uv,...
 for (int i = 0; i < numberOfParameterTypes; i++){
        allParameters[i]= std::unique_ptr<DataT[]>(new DataT[numberOfAllDataPoints*mFdimensions]);//To-Do:Fdim!!
        for( int j=0; j<numberOfAllDataPoints*mFdimensions; j++ ){
          allParameters[i][j] = -100;
        }
  }
  //filling allParameters[0] and FParameters with s:
  for (int i = 0; i < numberOfAllDataPoints; i++){
      allParameters[0][i]= DataPointF[i]; // TO DO - Just get the pointer adress there PLEASE!
      int p0indices[mXdimensions];
      arraytopoints(i,p0indices,numberOfDataPoints, mXdimensions);
      bool isKnot = 1;
        for (int j = 0; j < mXdimensions; j++){// is the current datapoint a knot?
          if (!mHelpers[j].getDataPoint(p0indices[j]).isKnot){isKnot = 0;break;} 
        }
        if (isKnot) { 
          int knotindices [mXdimensions];
          for (int j = 0; j < mXdimensions;j++ ){ // calculate KNotindices for all dimensions
            //WORKAROUND Getting Knotindices:
            knotindices[j]=p0indices[j]/((numberOfDataPoints[j]-1)/(numberOfKnots[j]-1));
            //knotindices[j] = mHelpers[j].getDataPoint(p0indices[j]).iKnot; //in der Annahme der wert ist ein Knotenindex und falls der datapoint ein knoten ist, gibt er seinen eigenen knotenindex zurück
          }
          // get the knotindexvalue for FParameters:
          int knotind = pointstoarray(knotindices, numberOfKnots, mXdimensions);

          Fparameters[knotind * numberOfParameterTypes] = DataPointF[i]; ///write derivatives in FParameters
        } // end if isKnot
  } //end i (filling DataPointF Values into allParameters[0] and FParameters)
  //now: allParameters[0] = dataPointF;
  
  //Array for input DataPointF-values for Spline1D::approximateFunctionGradually(...);
  std::unique_ptr<DataT[]> dataPointF1D[mXdimensions]; 
  for (int i = 0; i < mXdimensions; i++){
        dataPointF1D[i]= std::unique_ptr<DataT[]>(new DataT[numberOfDataPoints[i]]); // To-Do:Fdim!! For s and derivetives at all knots.
  }
  //Array to be filled by Spline1D::approximateFunctionGradually(...);
  std::unique_ptr<DataT[]> par[mXdimensions]; 
  for (int i = 0; i < mXdimensions; i++){
        par[i]= std::unique_ptr<DataT[]>(new DataT[numberOfKnots[i]*2]);
   }

  std::cout << "NumberOfParameters: " <<  mNumberOfParameters <<std::endl; 

  //STARTING MAIN-LOOP, for all Parametertypes:
  for (int p = 1; p < numberOfParameterTypes; p++ ){ // p = 1!! Wir kriegen s (p0) durch approximateFunction()oben
    int dimension; // find the dimension for approximation
    for (int i = (int)( log2f ((float) p)); i >= 0; i--){ 
      if (p%(int)(pow(2.0, i))==0){ dimension = i; break; }
    }
    dimension = ( p==1 || p==3 ) ?0 :1; //SG!!

    int currentDataPointF = p-(int)(pow(2.0, dimension));
    
    std::cout << std::endl << "p:" << p << ", dim of approximation: " << dimension << ", based on: " << currentDataPointF << std::endl;

    int nrOf1DSplines = (numberOfAllDataPoints / numberOfDataPoints[dimension]); // number of Splines for Parametertyp p in direction dim
    std::cout<<"N data points "<< numberOfAllDataPoints<<" N data points in dimension "<<numberOfDataPoints[dimension]<<std::endl;
    std::cout << "nr of splines: " << nrOf1DSplines;

    // getting the numbers of Datapoints for all dimension eccept the dimension of interpolation
    int currentNumbers[mXdimensions-1];
    for (int i = 0; i <dimension; i++){
      currentNumbers[i] = numberOfDataPoints[i];
    }   
    for (int i = dimension; i < mXdimensions-1; i++){
      currentNumbers[i] = numberOfDataPoints[i+1]; 
    }
    std::cout << " current numbers: ";
    for (int i = 0; i < mXdimensions-1; i++){
      std::cout << currentNumbers[i] << ",";
    }
    std::cout << std::endl;

    //// for all Splines in current dimension:
    for (int s = 0; s < nrOf1DSplines; s++){ 
      int indices[mXdimensions];
      arraytopoints  (s, indices, currentNumbers, mXdimensions-1);
      int startpoint [mXdimensions]; //startpoint for the current 1DSpline
      for (int i = 0; i<dimension; i++){
        startpoint[i]=indices[i];
      }
      startpoint[dimension] = 0;
      for (int i = dimension+1; i < mXdimensions; i++){
        startpoint[i] = indices[i-1];
      }
      
      bool isKnotOtherDim = 1;
      int knotindices [mXdimensions];
      for (int j = 0; j < mXdimensions; j++){//is dataPoint a knot?
        if( j== dimension ) continue;
        if (!mHelpers[j].getDataPoint(startpoint[j]).isKnot){isKnotOtherDim = 0;break;}
        //knotindices[j] = startpoint[j] / ((numberOfDataPoints[j]-1)/(numberOfKnots[j]-1));
        knotindices[j] = mHelpers[j].getDataPoint(startpoint[j]).iKnot; //in der Annahme der wert ist ein Knotenindex und falls der datapoint ein knoten ist, gibt er seinen eigenen knotenindex zurück
        if( startpoint[j] == numberOfDataPoints[j]-1 ) knotindices[j] = numberOfKnots[j]-1;
      }

      // NOW WE HAVE THE DATAPOINTINDICES OF THE CURRENT STARTPOINT IN startpoint-Array.
      int startdatapoint = pointstoarray( startpoint, numberOfDataPoints, mXdimensions);
      int distance = 1; // distance to the next dataPoint in the array for the current dimension
      for (int i = 0; i < dimension; i++){
        distance *= numberOfDataPoints[i]; 
      }
      distance *= 1; 

      for (int i = 0; i < numberOfDataPoints[dimension]; i++){ //Fill the dataPointF1D-Array
        dataPointF1D[dimension][i]   =    allParameters[currentDataPointF][startdatapoint + (i * distance)]; // uiuiui index kuddelmuddel???!!
      }
      mHelpers[dimension].approximateFunctionGradually(par[dimension].get(), dataPointF1D[dimension].get()); 
      // now we have all s and s' values in par[dimension]
      

      int redistributionindex[mXdimensions];
      for (int i = 0; i < mXdimensions; i++){
        redistributionindex [i] = startpoint[i];
      }
      //redistributing the derivatives at dimension-Knots into array p
      for (int i = 0; i < numberOfKnots[dimension];i++){ //for all dimension-Knots
        redistributionindex[dimension] = mHelpers[dimension].getKnotDataPoint(i);//find the indices
        int finalposition = pointstoarray(redistributionindex, numberOfDataPoints, mXdimensions);    
        allParameters[p][finalposition] = par[dimension][2*i +1];         
        if( isKnotOtherDim ){
          knotindices[dimension]=i;
          // get the knotindexvalue for FParameters:
          int knotind = pointstoarray(knotindices, numberOfKnots, mXdimensions);
          Fparameters[knotind * numberOfParameterTypes + p] = par[dimension][2*i +1]; ///write derivatives in FParameters
          std::cout<<"ND: parameter "<<p<<" ("<<redistributionindex[0]<<","<<redistributionindex[1]<<") = " << par[dimension][2*i +1]
           <<std::endl;
        } 
      } // end for all fknots (for redistribution) 

      // recalculation:
      
       for (int i = 0; i < numberOfDataPoints[dimension]; i++){// this is somehow still redundant// TO DO: ONLY PART OF approximateFunction WHERE NDIM is considerd!!
        redistributionindex[dimension] = i; // getting current datapointindices        
        if( !isKnotOtherDim ) continue; //SG!!!
        DataT splineF[1];
        DataT u = mHelpers[dimension].getDataPoint(i).u;
        mHelpers[dimension].getSpline().interpolateU(1, par[dimension].get(), u, splineF); //recalculate at all datapoints of dimension

        std::cout<<"recalculate point ("<<redistributionindex[0]<<","<<redistributionindex[1]<<")"
        <<" u "<<u<<" : "<<allParameters [currentDataPointF] [(int)(startdatapoint + i*distance )]
        <<" -> "<<splineF[0]
        <<std::endl;
        allParameters [currentDataPointF] [startdatapoint + i*distance] = splineF[0]; //write it in the array.                         
      } // end recalculation
      
    }//end of all1DSplines
  }//end of for parametertypes
} //end of approxymateFunction MYVERSION!  
 
  /*
  const int Ndim = mFdimensions;
  const int Ndim2 = 2 * Ndim;
  const int Ndim3 = 3 * Ndim;
  const int Ndim4 = 4 * Ndim;

  int nDataPointsU = getNumberOfDataPointsU1();
  int nDataPointsV = getNumberOfDataPointsU2();

  int nKnotsU = mHelperU1.getSpline().getNumberOfKnots();
  int nKnotsV = mHelperU2.getSpline().getNumberOfKnots();

  std::unique_ptr<DataT[]> rotDataPointF(new DataT[nDataPointsU * nDataPointsV * Ndim]); // U DataPoints x V DataPoints :  rotated DataPointF for one output dimension
  std::unique_ptr<DataT[]> Dv(new DataT[nKnotsV * nDataPointsU * Ndim]);                 // V knots x U DataPoints

  std::unique_ptr<DataT[]> parU(new DataT[mHelperU1.getSpline().getNumberOfParameters(Ndim)]);
  std::unique_ptr<DataT[]> parV(new DataT[mHelperU2.getSpline().getNumberOfParameters(Ndim)]);

  // rotated data points (u,v)->(v,u)

  for (int ipu = 0; ipu < nDataPointsU; ipu++) {
    for (int ipv = 0; ipv < nDataPointsV; ipv++) {
      for (int dim = 0; dim < Ndim; dim++) {
        rotDataPointF[Ndim * (ipu * nDataPointsV + ipv) + dim] = DataPointF[Ndim * (ipv * nDataPointsU + ipu) + dim];
      }
    }
  }

  // get S and S'u at all the knots by interpolating along the U axis

  for (int iKnotV = 0; iKnotV < nKnotsV; ++iKnotV) {
    int ipv = mHelperU2.getKnotDataPoint(iKnotV);
    const DataT* DataPointFrow = &(DataPointF[Ndim * ipv * nDataPointsU]);
    mHelperU1.approximateFunctionGradually(parU.get(), DataPointFrow);

    for (int iKnotU = 0; iKnotU < nKnotsU; ++iKnotU) {
      DataT* knotPar = &Fparameters[Ndim4 * (iKnotV * nKnotsU + iKnotU)];
      for (int dim = 0; dim < Ndim; ++dim) {
        knotPar[dim] = parU[Ndim * (2 * iKnotU) + dim];                // store S for all the knots
        knotPar[Ndim2 + dim] = parU[Ndim * (2 * iKnotU) + Ndim + dim]; // store S'u for all the knots //SG!!!
      }
    }

    // recalculate F values for all ipu DataPoints at V = ipv
    for (int ipu = 0; ipu < nDataPointsU; ipu++) {
      DataT splineF[Ndim];
      DataT u = mHelperU1.getDataPoint(ipu).u;
      mHelperU1.getSpline().interpolateU(Ndim, parU.get(), u, splineF);
      for (int dim = 0; dim < Ndim; dim++) {
        rotDataPointF[(ipu * nDataPointsV + ipv) * Ndim + dim] = splineF[dim];
      }
    }
  }

  // calculate S'v at all data points with V == V of a knot

  for (int ipu = 0; ipu < nDataPointsU; ipu++) {
    const DataT* DataPointFcol = &(rotDataPointF[ipu * nDataPointsV * Ndim]);
    mHelperU2.approximateFunctionGradually(parV.get(), DataPointFcol);
    for (int iKnotV = 0; iKnotV < nKnotsV; iKnotV++) {
      for (int dim = 0; dim < Ndim; dim++) {
        DataT dv = parV[(iKnotV * 2 + 1) * Ndim + dim];
        Dv[(iKnotV * nDataPointsU + ipu) * Ndim + dim] = dv;
      }
    }
  }

  // fit S'v and S''_vu at all the knots

  for (int iKnotV = 0; iKnotV < nKnotsV; ++iKnotV) {
    const DataT* Dvrow = &(Dv[iKnotV * nDataPointsU * Ndim]);
    mHelperU1.approximateFunction(parU.get(), Dvrow);
    for (int iKnotU = 0; iKnotU < nKnotsU; ++iKnotU) {
      for (int dim = 0; dim < Ndim; ++dim) {
        Fparameters[Ndim4 * (iKnotV * nKnotsU + iKnotU) + Ndim + dim] = parU[Ndim * 2 * iKnotU + dim];         // store S'v for all the knots
        Fparameters[Ndim4 * (iKnotV * nKnotsU + iKnotU) + Ndim3 + dim] = parU[Ndim * 2 * iKnotU + Ndim + dim]; // store S''vu for all the knots
      }
    }
  }
  */


template class GPUCA_NAMESPACE::gpu::SplineHelperTest<float>;
template class GPUCA_NAMESPACE::gpu::SplineHelperTest<double>;

#endif
