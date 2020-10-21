// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  SplineSpec.h
/// \brief Definition of SplineSpec class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#ifndef ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SPLINESPEC_H
#define ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SPLINESPEC_H

#include "Spline1D.h"
#include "FlatObject.h"
#include "GPUCommonDef.h"

#if !defined(__CINT__) && !defined(__ROOTCINT__) && !defined(GPUCA_GPUCODE) && !defined(GPUCA_NO_VC) && defined(__cplusplus) && __cplusplus >= 201703L
#include <Vc/Vc>
#include <Vc/SimdArray>
#endif

class TFile;

namespace GPUCA_NAMESPACE
{
namespace gpu
{

/// ==================================================================================================
/// The class SplineContainer is a base class of Spline.
/// It contains all the class members and methods which only depends on the DataT data type
/// and do not depend on other template parameters of Spline.
/// It also contains all non-inlined methods with the implementation in SplineSpec.cxx file.
///
/// DataT is a data type, which is supposed to be either double or float.
/// For other possible data types one has to add the corresponding instantiation line
/// at the end of the SplineSpec.cxx file
///
template <typename DataT>
class SplineContainer : public FlatObject
{
 public:
  typedef typename Spline1D<DataT>::SafetyLevel SafetyLevel;
  typedef typename Spline1D<DataT>::Knot Knot;

  /// _____________  Version control __________________________

  /// Version control
  GPUd() static constexpr int getVersion() { return (1 << 16) + Spline1D<DataT>::getVersion(); }

  /// _____________  C++ constructors / destructors __________________________

  /// Default constructor
  SplineContainer() CON_DEFAULT;

  /// Disable all other constructors
  SplineContainer(const SplineContainer&) CON_DELETE;

  /// Destructor
  ~SplineContainer() CON_DEFAULT;

  /// _______________  Construction interface  ________________________

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE)
  /// approximate a function F with this spline
  void approximateFunction(const double xMin[/* mXdim */], const double xMax[/* mXdim */],
                           std::function<void(const double x[/* mXdim */], double f[/*mYdim*/])> F,
                           const int nAuxiliaryDataPoints[/* mXdim */] = nullptr);
#endif

  /// _______________  IO   ________________________

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE)
  /// write a class object to the file
  int writeToFile(TFile& outf, const char* name);

  /// read a class object from the file
  static SplineContainer* readFromFile(TFile& inpf, const char* name);
#endif

  /// _______________  Getters   ________________________

  /// Get number of X dimensions
  GPUd() int getXdimensions() const { return mXdim; }

  /// Get number of Y dimensions
  GPUd() int getYdimensions() const { return mYdim; }

  /// Get minimal required alignment for the spline parameters
  GPUd() static constexpr size_t getParameterAlignmentBytes() { return 16; }

  /// Number of parameters
  GPUd() int getNumberOfParameters() const { return this->calcNumberOfParameters(mYdim); }

  /// Size of the parameter array in bytes
  GPUd() size_t getSizeOfParameters() const { return sizeof(DataT) * this->getNumberOfParameters(); }

  /// Get a number of knots
  GPUd() int getNumberOfKnots() const { return mNknots; }

  /// Number of parameters per knot
  GPUd() int getNumberOfParametersPerKnot() const { return calcNumberOfParametersPerKnot(mYdim); }

  /// Get 1-D grid for dimX dimension
  GPUd() const Spline1D<DataT>& getGrid(int dimX) const { return mGrid[dimX]; }

  /// Get u[] coordinate of i-th knot
  GPUd() void getKnotU(int iKnot, int u[/* mXdim */]) const;

  /// Get index of a knot (iKnot1,iKnot2,..,iKnotN)
  GPUd() int getKnotIndex(const int iKnot[/* mXdim */]) const;

  /// Get spline parameters
  GPUd() DataT* getParameters() { return mParameters; }

  /// Get spline parameters const
  GPUd() const DataT* getParameters() const { return mParameters; }

  /// _______________  Technical stuff  ________________________

  /// Get offset of Grid[dimX] flat data in the flat buffer
  GPUd() size_t getGridOffset(int dimX) const { return mGrid[dimX].getFlatBufferPtr() - mFlatBufferPtr; }

  /// Set X range
  void setXrange(const DataT xMin[/* mXdim */], const DataT xMax[/* mXdim */]);

  /// Print method
  void print() const;

  ///  _______________  Expert tools  _______________

  /// Number of parameters for given Y dimensions
  GPUd() int calcNumberOfParameters(int nYdim) const
  {
    return calcNumberOfParametersPerKnot(nYdim) * getNumberOfKnots();
  }

  /// Number of parameters per knot
  GPUd() int calcNumberOfParametersPerKnot(int nYdim) const
  {
    return (1 << mXdim) * nYdim; // 2^mXdim parameters per Y dimension
  }

    ///_______________  Test tools  _______________

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE) // code invisible on GPU and in the standalone compilation
  /// Test the class functionality
  static int test(const bool draw = 0, const bool drawDataPoints = 1);
#endif

  /// _____________  FlatObject functionality, see FlatObject class for description  ____________

  using FlatObject::getBufferAlignmentBytes;
  using FlatObject::getClassAlignmentBytes;

#if !defined(GPUCA_GPUCODE)
  void cloneFromObject(const SplineContainer& obj, char* newFlatBufferPtr);
  void moveBufferTo(char* newBufferPtr);
#endif

  using FlatObject::releaseInternalBuffer;

  void destroy();
  void setActualBufferAddress(char* actualFlatBufferPtr);
  void setFutureBufferAddress(char* futureFlatBufferPtr);

 protected:
#if !defined(GPUCA_GPUCODE)
  /// Constructor for a regular spline
  void recreate(int nXdim, int nYdim, const int nKnots[/* nXdim */]);

  /// Constructor for an irregular spline
  void recreate(int nXdim, int nYdim, const int nKnots[/* nXdim */], const int* const knotU[/* nXdim */]);
#endif

  /// _____________  Data members  ____________

  int mXdim = 0;   ///< dimentionality of X
  int mYdim = 0;   ///< dimentionality of Y
  int mNknots = 0; ///< number of spline knots

  Spline1D<DataT>* mGrid; //! (transient!!) mXdim grids
  DataT* mParameters;     //! (transient!!) F-dependent parameters of the spline

#ifndef GPUCA_ALIROOT_LIB
  ClassDefNV(SplineContainer, 1);
#endif
};

template <typename DataT>
GPUdi() void SplineContainer<DataT>::getKnotU(int iKnot, int u[/* mXdim */]) const
{
  /// Get u[] coordinate of i-th knot
  for (int dim = 0; dim < mXdim; dim++) {
    int n = mGrid[dim].getNumberOfKnots();
    u[dim] = mGrid[dim].getKnot(iKnot % n).getU();
    iKnot /= n;
  }
}

template <typename DataT>
GPUdi() int SplineContainer<DataT>::getKnotIndex(const int iKnot[/* mXdim */]) const
{
  /// Get index of a knot (iKnot1,iKnot2,..,iKnotN)
  int ind = iKnot[0];
  int n = 1;
  for (int dim = 1; dim < mXdim; dim++) {
    n *= mGrid[dim - 1].getNumberOfKnots();
    ind += n * iKnot[dim];
  }
  return ind;
}

template <typename DataT>
GPUdi() void SplineContainer<DataT>::
  setXrange(const DataT xMin[/* mXdim */], const DataT xMax[/* mXdim */])
{
  /// Set X range
  for (int i = 0; i < mXdim; i++) {
    mGrid[i].setXrange(xMin[i], xMax[i]);
  }
}

/// ==================================================================================================
/// The specification declares common methods for all other Spline specifications.
/// Implementations of the methods may depend on the YdimT value.
///
template <typename DataT, int XdimT, int YdimT>
class SplineCommon : public SplineContainer<DataT>
{
  typedef SplineContainer<DataT> TBase;

 public:
  typedef typename TBase::SafetyLevel SafetyLevel;
  typedef typename TBase::Knot Knot;

  /// _______________  Template magic  ________________________

  /// An expression getXdim(XdimCase{},int) is either an integer or a constexpr integer,
  /// depending on the (XdimT>0) value
  ///
  typedef std::integral_constant<bool, (XdimT > 0)> XdimCase;
  static constexpr int getXdim(std::true_type, int) { return XdimT; }
  static int getXdim(std::false_type, int nXdim) { return nXdim; }

  typedef std::integral_constant<bool, (YdimT > 0)> YdimCase;
  static constexpr int getYdim(std::true_type, int) { return YdimT; }
  static int getYdim(std::false_type, int nYdim) { return nYdim; }

  /// An expression getMaxXdim(XdimMaxCase{},int) is either an integer or a constexpr integer,
  /// depending on the (XdimT==0) value
  ///
  typedef std::integral_constant<bool, (XdimT == 0)> XdimMaxCase;
  static int getMaxXdim(std::true_type, int nXdim) { return nXdim; }
  static constexpr int getMaxXdim(std::false_type, int) { return abs(XdimT); }

  typedef std::integral_constant<bool, (YdimT == 0)> YdimMaxCase;
  static int getMaxYdim(std::true_type, int nYdim) { return nYdim; }
  static constexpr int getMaxYdim(std::false_type, int) { return abs(YdimT); }

  /// _______________  Interpolation math   ________________________

  /// Get interpolated value S(x)
  GPUd() void interpolate(const DataT x[/*mXdim*/], GPUgeneric() DataT S[/*mYdim*/]) const
  {
    auto nXdim = getXdim(XdimCase{}, mXdim);
    auto maxXdim = getMaxXdim(XdimMaxCase{}, mXdim);
    DataT u[maxXdim];
    for (int i = 0; i < nXdim; i++) {
      u[i] = mGrid[i].convXtoU(x[i]);
    }
    interpolateU<SafetyLevel::kSafe>(mXdim, mYdim, mParameters, u, S);
  }

  /// Get interpolated value for S(u):inpXdim->inpYdim using spline parameters Parameters
  template <SafetyLevel SafeT = SafetyLevel::kSafe>
  GPUd() void interpolateU(int inpXdim, int inpYdim, GPUgeneric() const DataT Parameters[],
                           const DataT u[/*inpXdim*/], GPUgeneric() DataT S[/*inpYdim*/]) const
  {
    auto nXdim = getXdim(XdimCase{}, inpXdim);
    auto nYdim = getYdim(YdimCase{}, inpYdim);
    auto maxXdim = getMaxXdim(XdimMaxCase{}, inpXdim);
    auto maxYdim = getMaxYdim(YdimMaxCase{}, inpYdim);

    auto nParameters = 1 << (2 * nXdim);         //total Nr of Parameters necessary for one interpolation
    auto nKnotParametersPerY = 1 << nXdim;       // Nr of Parameters per Knot per Y dimension
    auto nKnotParameters = (1 << nXdim) * nYdim; // Nr of Parameters per Knot

    DataT iParameters[(1 << (2 * maxXdim)) * maxYdim]; // Array for all parameters

    //get the indices of the "most left" Knot:

    int indices[maxXdim]; //indices of the 'most left' knot
    for (int i = 0; i < nXdim; i++) {
      indices[i] = mGrid[i].getLeftKnotIndexForU(u[i]);
    }
    // get all the needed parameters into one array iParameters[nParameters]:
    int indicestmp[maxXdim];
    for (int i = 0; i < nKnotParametersPerY; i++) { // for every necessary Knot
      for (int k = 0; k < nXdim; k++) {
        indicestmp[k] = indices[k] + (i / (1 << k)) % 2; //get the knot-indices in every dimension (mirrored order binary counting)
      }
      int index = TBase::getKnotIndex(indicestmp); //get index of the current Knot

      for (int j = 0; j < nKnotParameters; j++) { //and fill the iparameter array with according parameters
        iParameters[i * nKnotParameters + j] = Parameters[index * nKnotParameters + j];
      }
    }
    //now start with the interpolation loop:

    auto maxInterpolations = (1 << (2 * maxXdim - 2)) * maxYdim;

    DataT S0[maxInterpolations];
    DataT D0[maxInterpolations];
    DataT S1[maxInterpolations];
    DataT D1[maxInterpolations];

    int nrofInterpolations = (1 << (2 * nXdim - 2)) * nYdim;
    int nrofKnots = 1 << (nXdim);

    for (int d = 0; d < nXdim; d++) {            //for every dimension
      DataT* pointer[4] = {S0, D0, S1, D1};      // pointers for interpolation arrays S0, D0, S1, D1 point to Arraystart
      for (int i = 0; i < nrofKnots; i++) {      //for every knot
        for (int j = 0; j < nrofKnots; j++) {    // for every parametertype
          int pointernr = 2 * (i % 2) + (j % 2); //to which array should it be delivered
          for (int k = 0; k < nYdim; k++) {
            pointer[pointernr][0] = iParameters[(i * nrofKnots + j) * nYdim + k];
            pointer[pointernr]++;
          }
        } // end for j (every parametertype)
      }   // end for i (every knot)

      const typename Spline1D<DataT>::Knot& knotL = mGrid[d].getKnot(indices[d]);
      int Ydim = nrofInterpolations;
      DataT coordinate = u[d];

      typedef Spline1DSpec<DataT, YdimT, true, (YdimT > 0), (YdimT == 1), (YdimT == 999999)> TGridX;
      const TGridX& gridX = *((const TGridX*)&(mGrid[d]));
      gridX.interpolateU(Ydim, knotL, S0, D0, S1, D1, coordinate, iParameters);

      nrofInterpolations = nrofInterpolations / 4;
      nrofKnots = nrofKnots / 2;
    } //end d (every dimension)

    for (int i = 0; i < nYdim; i++) {
      S[i] = iParameters[i]; // write into result-array
      //std::cout<<iParameters[i] <<", ";
    }
  } // end interpolateU

 protected:
  using TBase::TBase; // inherit constructors and hide them
  using TBase::mGrid;
  using TBase::mParameters;
  using TBase::mXdim;
  using TBase::mYdim;
};

/// SplineSpec class declares different specifications of the Spline class.
/// (See Spline.h for the description.)
///
/// The specifications depend on the value of Spline's template parameters XdimT and YdimT.
/// Specifications have different constructors and slightly different declarations of methods.
///
/// The meaning of the template parameters:
///
/// \param DataT data type: float or double
/// \param XdimT
///    XdimT > 0 : the number of X dimensions is known at the compile time and is equal to XdimT
///    XdimT = 0 : the number of X dimensions will be set in the runtime
///    XdimT < 0 : the number of X dimensions will be set in the runtime, and it will not exceed abs(XdimT)
/// \param XdimTPositive A technical flag == (XdimT > 0). It lets us to have different specializations
///                      for XdimT>0 and XdimT<=0 regions
/// \param XdimTGeneric  A technical flag. It lets an explicit specialization (for a particular XdimT value)
///                          to inherit from the generic specialization (for any XdimT)
/// \param YdimT same meaning as YdimT
/// \param YdimTPositive same meaning as XdimTPositive
/// \param YdimTGeneric same meaning as XdimTMarkedGeneric
///
template <typename DataT,
          int XdimT, bool XdimTPositive, bool XdimTGeneric,
          int YdimT, bool YdimTPositive, bool YdimTGeneric>
class SplineSpec;

/// ==================================================================================================
/// A specification (XdimT<=0 || YdimT<=0) where at least one of the dimensions
/// must be set in the runtime via a constructor parameter
///
template <typename DataT,
          int XdimT, bool XdimTPositive, bool XdimTGeneric,
          int YdimT, bool YdimTPositive, bool YdimTGeneric>
class SplineSpec
  : public SplineCommon<DataT, XdimT, YdimT>
{
  typedef SplineContainer<DataT> TVeryBase;
  typedef SplineCommon<DataT, XdimT, YdimT> TBase;

 public:
  typedef typename TVeryBase::SafetyLevel SafetyLevel;

#if !defined(GPUCA_GPUCODE)
  /// Default constructor
  SplineSpec() : SplineSpec((XdimT > 0 ? XdimT : 0), (YdimT > 0 ? YdimT : 0), nullptr) {}

  /// Constructor for a regular spline
  SplineSpec(int nXdim, int nYdim, const int nKnots[/* nXdim */]) : TBase()
  {
    this->recreate(nXdim, nYdim, nKnots);
  }

  /// Constructor for an irregular spline
  SplineSpec(int nXdim, int nYdim, const int nKnots[/* nXdim */], const int* const knotU[/* nXdim */])
    : TBase()
  {
    this->recreate(nXdim, nYdim, nKnots, knotU);
  }

  /// Copy constructor
  SplineSpec(const SplineSpec& v) : TBase()
  {
    cloneFromObject(v, nullptr);
  }

  /// Constructor for a regular spline
  void recreate(int nXdim, int nYdim, const int nKnots[/* nXdim */])
  {
    checkDimensions(nXdim, nYdim);
    TBase::recreate(nXdim, nYdim, nKnots);
  }

  /// Constructor for an irregular spline
  void recreate(int nXdim, int nYdim, const int nKnots[/* nXdim */], const int* const knotU[/* nXdim */])
  {
    checkDimensions(nXdim, nYdim);
    TBase::recreate(nXdim, nYdim, nKnots, knotU);
  }

#endif

  ///  _______  Expert tools: interpolation with given nYdim and external Parameters _______

  using TBase::interpolateU;

  /// Check dimensions
  void checkDimensions(int& nXdim, int& nYdim)
  {
    if (XdimT > 0 && nXdim != XdimT) {
      assert(0);
      nXdim = XdimT;
    }
    if (XdimT < 0 && nXdim > abs(XdimT)) {
      assert(0);
      nXdim = abs(XdimT);
    }
    if (nXdim < 0) {
      assert(0);
      nXdim = 0;
    }
    if (YdimT > 0 && nYdim != YdimT) {
      assert(0);
      nYdim = YdimT;
    }
    if (YdimT < 0 && nYdim > abs(YdimT)) {
      assert(0);
      nYdim = abs(YdimT);
    }
    if (nYdim < 0) {
      assert(0);
      nYdim = 0;
    }
  }
};

/// ==================================================================================================
/// Specification XdimT>0, YdimT>0 where the number of dimensions is taken from template parameters
/// at the compile time
///
template <typename DataT,
          int XdimT, bool XdimTGeneric,
          int YdimT, bool YdimTGeneric>
class SplineSpec<DataT,
                 XdimT, true /*XdimTPositive*/, XdimTGeneric,
                 YdimT, true /*YdimTPositive*/, YdimTGeneric>
  : public SplineCommon<DataT, XdimT, YdimT>
{
  typedef SplineContainer<DataT> TVeryBase;
  typedef SplineCommon<DataT, XdimT, YdimT> TBase;

 public:
  typedef typename TVeryBase::SafetyLevel SafetyLevel;

#if !defined(GPUCA_GPUCODE)
  /// Default constructor
  SplineSpec() : SplineSpec(nullptr) {}

  /// Constructor for a regular spline
  SplineSpec(const int nKnots[/*XdimT*/]) : TBase()
  {
    recreate(nKnots);
  }
  /// Constructor for an irregular spline
  SplineSpec(const int nKnots[/*XdimT*/], const int* const knotU[/*XdimT*/])
    : TBase()
  {
    recreate(nKnots, knotU);
  }
  /// Copy constructor
  SplineSpec(const SplineSpec& v) : TBase()
  {
    TBase::cloneFromObject(v, nullptr);
  }
  /// Constructor for a regular spline
  void recreate(const int nKnots[/*XdimT*/])
  {
    TBase::recreate(XdimT, YdimT, nKnots);
  }

  /// Constructor for an irregular spline
  void recreate(const int nKnots[/*XdimT*/], const int* const knotU[/*XdimT*/])
  {
    TBase::recreate(XdimT, YdimT, nKnots, knotU);
  }
#endif

  /// Get number of X dimensions
  GPUd() constexpr int getXdimensions() const { return XdimT; }

  /// Get number of Y dimensions
  GPUd() constexpr int getYdimensions() const { return YdimT; }

  ///  _______  Expert tools: interpolation with given nYdim and external Parameters _______

  /// Get interpolated value for an YdimT-dimensional S(u1,u2) using spline parameters Parameters.
  template <SafetyLevel SafeT = SafetyLevel::kSafe>
  GPUd() void interpolateU(GPUgeneric() const DataT Parameters[],
                           const DataT u[/*XdimT*/], GPUgeneric() DataT S[/*YdimT*/]) const
  {
    TBase::template interpolateU<SafeT>(XdimT, YdimT, Parameters, u, S);
  }

  /// _______________  Suppress some parent class methods   ________________________
 private:
#if !defined(GPUCA_GPUCODE)
  using TBase::recreate;
#endif
  using TBase::interpolateU;
};

/// ==================================================================================================
/// Specification where the number of Y dimensions is 1.
///
template <typename DataT,
          int XdimT, bool XdimTPositive, bool XdimTGeneric>
class SplineSpec<DataT,
                 XdimT, XdimTPositive, XdimTGeneric,
                 1, true /*YdimTPositive*/, false /*YdimTGeneric*/>
                 {};
#ifdef XXX
  : public SplineSpec<DataT,
                      XdimT, XdimTPositive, XdimTGeneric,
                      1, true /*YdimTPositive*/, true /*YdimTGeneric*/>
{
  typedef SplineSpec<DataT,
                     XdimT, XdimTPositive, XdimTGeneric,
                     1, true /*YdimTPositive*/, true /*YdimTGeneric*/>
    TBase;

 public:
  using TBase::TBase; // inherit constructors

  /// Simplified interface for 1D: return the interpolated value
  GPUd() DataT interpolate(const DataT x[]) const
  {
    DataT S = 0.;
    TBase::interpolate(x, &S);
    return S;
  }

  // this parent method should be public anyhow,
  // but the compiler gets confused w/o this extra declaration
  using TBase::interpolate;
};
#endif

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
