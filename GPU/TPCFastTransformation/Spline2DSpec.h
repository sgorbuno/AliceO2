// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  Spline2DSpec.h
/// \brief Definition of Spline2DSpec class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#ifndef ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SPLINE2DSPEC_H
#define ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SPLINE2DSPEC_H

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
/*
/// Spline2DSpec class declares different specifications of the Spline2D class.
/// (See Spline2D.h for the description.)
///
/// The specifications depend on the value of Spline2D's template parameter YdimT.
/// Specifications have different constructors and slightly different declarations of methods.
///
/// The meaning of the template parameters:
///
/// \param DataT data type: float or double
/// \param YdimT
///     >= 0 : YdimT is the number of Y dimensions. Use it when it is known at the compile time.
///  not set : the number of Y dimensions not known and will be set in the runtime
///     < 0  : the number of Y dimensions will be set in the runtime, but it will not exceed abs(YdimT)
/// \param YisAnyT      YdimT is any. This case is a parent for all other specifications.
/// \param YisPositiveT YdimT >= 0
/// \param YisOneT      YdimT == 1
/// \param YisNotSetT   YdimT is not set (it is equal to some specific default value)
///
template <typename DataT, int YdimT, bool YisAnyT, bool YisPositiveT, bool YisOneT, bool YisNotSetT>
class Spline2DSpec;
*/

/// ==================================================================================================
/// The class Spline2DContainer is a base class of Spline2D.
/// It contains all the class members and methods which only depends on the DataT data type
/// and do not depend on other template parameters of Spline2D.
/// It also contains all non-inlined methods with the implementation in Spline2DSpec.cxx file.
///
/// DataT is a data type, which is supposed to be either double or float.
/// For other possible data types one has to add the corresponding instantiation line
/// at the end of the Spline2DSpec.cxx file
///
template <typename DataT>
class Spline2DContainer : public FlatObject
{
 public:
  typedef typename Spline1D<DataT>::SafetyLevel SafetyLevel;
  typedef typename Spline1D<DataT>::Knot Knot;

  /// _____________  Version control __________________________

  /// Version control
  GPUd() static constexpr int getVersion() { return (1 << 16) + Spline1D<DataT>::getVersion(); }

  /// _____________  C++ constructors / destructors __________________________

  /// Default constructor
  Spline2DContainer() CON_DEFAULT;

  /// Disable all other constructors
  Spline2DContainer(const Spline2DContainer&) CON_DELETE;

  /// Destructor
  ~Spline2DContainer() CON_DEFAULT;

/// _______________  Construction interface  ________________________

#if !defined(GPUCA_GPUCODE)
  /// approximate a function F with this spline
  void approximateFunction(double x1Min, double x1Max, double x2Min, double x2Max,
                           std::function<void(double x1, double x2, double f[/*mYdim*/])> F,
                           int nAuxiliaryDataPointsU1 = 4, int nAuxiliaryDataPointsU2 = 4);
#endif

/// _______________  IO   ________________________

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE)
/*
  /// write a class object to the file
  int writeToFile(TFile& outf, const char* name);

  /// read a class object from the file
  static Spline2DContainer* readFromFile(TFile& inpf, const char* name);
  */
#endif

  /// _______________  Getters   ________________________

  /// Get number of Y dimensions
  GPUd() int getYdimensions() const { return mYdim; }

  /// Get minimal required alignment for the spline parameters
  GPUd() static constexpr size_t getParameterAlignmentBytes() { return 16; }

  /// Number of parameters
  GPUd() int getNumberOfParameters() const { return this->calcNumberOfParameters(mYdim); }

  /// Size of the parameter array in bytes
  GPUd() size_t getSizeOfParameters() const { return sizeof(DataT) * this->getNumberOfParameters(); }

  /// Get a number of knots
  GPUd() int getNumberOfKnots() const { return mGridX1.getNumberOfKnots() * mGridX2.getNumberOfKnots(); }

  /// Get 1-D grid for the X1 coordinate
  GPUd() const Spline1D<DataT>& getGridX1() const { return mGridX1; }

  /// Get 1-D grid for the X2 coordinate
  GPUd() const Spline1D<DataT>& getGridX2() const { return mGridX2; }

  /// Get 1-D grid for X1 or X2 coordinate
  GPUd() const Spline1D<DataT>& getGrid(int ix) const { return (ix == 0) ? mGridX1 : mGridX2; }

  /// Get (u1,u2) of i-th knot
  GPUd() void getKnotU(int iKnot, int& u1, int& u2) const
  {
    u1 = mGridX1.getKnot(iKnot % mGridX1.getNumberOfKnots()).getU();
    u2 = mGridX2.getKnot(iKnot / mGridX1.getNumberOfKnots()).getU();
  }

  /// Get index of a knot (iKnotX1,iKnotX2)
  GPUd() int getKnotIndex(int iKnotX1, int iKnotX2) const
  {
    return mGridX1.getNumberOfKnots() * iKnotX2 + iKnotX1;
  }

  /// Get spline parameters
  GPUd() DataT* getParameters() { return mParameters; }

  /// Get spline parameters const
  GPUd() const DataT* getParameters() const { return mParameters; }

  /// _______________  Technical stuff  ________________________

  /// Get offset of GridX1 flat data in the flat buffer
  GPUd() size_t getGridX1Offset() const { return mGridX1.getFlatBufferPtr() - mFlatBufferPtr; }

  /// Get offset of GridX2 flat data in the flat buffer
  GPUd() size_t getGridX2Offset() const { return mGridX2.getFlatBufferPtr() - mFlatBufferPtr; }

  /// Set X range
  GPUd() void setXrange(DataT x1Min, DataT x1Max, DataT x2Min, DataT x2Max)
  {
    mGridX1.setXrange(x1Min, x1Max);
    mGridX2.setXrange(x2Min, x2Max);
  }

  /// Print method
  void print() const;

  ///  _______________  Expert tools  _______________

  /// Number of parameters for given Y dimensions
  GPUd() int calcNumberOfParameters(int nYdim) const { return (4 * nYdim) * getNumberOfKnots(); }

///_______________  Test tools  _______________

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE) // code invisible on GPU and in the standalone compilation
  /// Test the class functionality
  static int test(const bool draw = 0, const bool drawDataPoints = 1);
#endif

  /// _____________  FlatObject functionality, see FlatObject class for description  ____________

  using FlatObject::getBufferAlignmentBytes;
  using FlatObject::getClassAlignmentBytes;

#if !defined(GPUCA_GPUCODE)
  void cloneFromObject(const Spline2DContainer& obj, char* newFlatBufferPtr);
  void moveBufferTo(char* newBufferPtr);
#endif

  using FlatObject::releaseInternalBuffer;

  void destroy();
  void setActualBufferAddress(char* actualFlatBufferPtr);
  void setFutureBufferAddress(char* futureFlatBufferPtr);

 protected:
#if !defined(GPUCA_GPUCODE)
  /// Constructor for a regular spline
  void recreate(int nYdim, int numberOfKnotsU1, int numberOfKnotsU2);

  /// Constructor for an irregular spline
  void recreate(int nYdim, int numberOfKnotsU1, const int knotU1[], int numberOfKnotsU2, const int knotU2[]);
#endif

  /// _____________  Data members  ____________

  int mYdim = 0;                ///< dimentionality of F
  Spline1D<DataT> mGridX1;      ///< grid for U axis
  Spline1D<DataT> mGridX2;      ///< grid for V axis
  DataT* mParameters = nullptr; //! (transient!!) F-dependent parameters of the spline

#ifndef GPUCA_ALIROOT_LIB
  ClassDefNV(Spline2DContainer, 1);
#endif
};

#ifdef XXX
/// ==================================================================================================
/// Specification (YisAnyT==1, YisPositiveT=*, YisOneT=*, YisNotSetT=*)
/// It  declares common methods for all other Spline2D specifications.
/// Implementations may depend on the YdimT value.
///
template <typename DataT, int YdimT, bool YisPositiveT, bool YisOneT, bool YisNotSetT>
class Spline2DSpec<DataT, YdimT, true, YisPositiveT, YisOneT, YisNotSetT> : public Spline2DContainer<DataT>
{
  typedef Spline2DContainer<DataT> TBase;
  typedef typename TBase::SafetyLevel SafetyLevel;
  typedef typename TBase::Knot Knot;

 public:
  /// _______________  Main functionality   ________________________

  /// Get interpolated value S(x)
  GPUd() void interpolate(DataT x1, DataT x2, GPUgeneric() DataT S[/*mYdim*/]) const
  {
    interpolateU<SafetyLevel::kSafe>(mYdim, mParameters, mGridX1.convXtoU(x1), mGridX2.convXtoU(x2), S);
  }

 protected:
  using TBase::TBase; // inherit constructors and hide them
  using TBase::mGridX1;
  using TBase::mGridX2;
  using TBase::mParameters;
  using TBase::mYdim;

  /// A template magic.
  /// An expression getYdim(dimType{},int) is either an integer or a constexpr integer,
  /// depending on the YisPositiveT value
  ///
  static constexpr int getYdim(std::true_type, int) { return YdimT; }
  static int getYdim(std::false_type, int nYdim) { return nYdim; }
  typedef std::integral_constant<bool, YisPositiveT> dimType;

  /// Get interpolated value for an inpYdim-dimensional S(u1,u2) using spline parameters Parameters.
  template <SafetyLevel SafeT = SafetyLevel::kSafe>
  GPUd() void interpolateU(int inpYdim, GPUgeneric() const DataT Parameters[],
                           DataT u1, DataT u2, GPUgeneric() DataT S[/*inpYdim*/]) const
  {
    auto nYdim = getYdim(dimType{}, inpYdim);

    /*
  float u = u1;
  float v = u2;
  int nu = mGridU1.getNumberOfKnots();
  int iu = mGridU1.getLeftKnotIndexForU<Safe>(u);
  int iv = mGridU2.getLeftKnotIndexForU<Safe>(v);

  const typename Spline1D<DataT>::Knot& knotU = mGridU1.getKnot<SafetyLevel::kNotSafe>(iu);
  const typename Spline1D<DataT>::Knot& knotV = mGridU2.getKnot<SafetyLevel::kNotSafe>(iv);

  const int nYdim2 = nYdim * 2;
  const int nYdim4 = nYdim * 4;

  const DataT* par00 = Parameters + (nu * iv + iu) * nYdim4; // values { {Y1,Y2,Y3}, {Y1,Y2,Y3}'v, {Y1,Y2,Y3}'u, {Y1,Y2,Y3}''vu } at {u0, v0}
  const DataT* par10 = par00 + nYdim4;                        // values { ... } at {u1, v0}
  const DataT* par01 = par00 + nYdim4 * nu;                   // values { ... } at {u0, v1}
  const DataT* par11 = par01 + nYdim4;                        // values { ... } at {u1, v1}

  DataT Su0[nYdim4]; // values { {Y1,Y2,Y3,Y1'v,Y2'v,Y3'v}(v0), {Y1,Y2,Y3,Y1'v,Y2'v,Y3'v}(v1) }, at u0
  DataT Du0[nYdim4]; // derivatives {}'_u  at u0
  DataT Su1[nYdim4]; // values { {Y1,Y2,Y3,Y1'v,Y2'v,Y3'v}(v0), {Y1,Y2,Y3,Y1'v,Y2'v,Y3'v}(v1) }, at u1
  DataT Du1[nYdim4]; // derivatives {}'_u  at u1

  for (int i = 0; i < nYdim2; i++) {
    Su0[i] = par00[i];
    Su0[nYdim2 + i] = par01[i];

    Du0[i] = par00[nYdim2 + i];
    Du0[nYdim2 + i] = par01[nYdim2 + i];

    Su1[i] = par10[i];
    Su1[nYdim2 + i] = par11[i];

    Du1[i] = par10[nYdim2 + i];
    Du1[nYdim2 + i] = par11[nYdim2 + i];
  }

  DataT parU[nYdim4]; // interpolated values { {Y1,Y2,Y3,Y1'v,Y2'v,Y3'v}(v0), {Y1,Y2,Y3,Y1'v,Y2'v,Y3'v}(v1) } at u
  mGridU1.interpolateUMath(nYdim4, knotU, Su0, Du0, Su1, Du1, u, parU);

  const DataT* Sv0 = parU + 0;
  const DataT* Dv0 = parU + nYdim;
  const DataT* Sv1 = parU + nYdim2;
  const DataT* Dv1 = parU + nYdim2 + nYdim;
  if (DoFirstDimOnly) {
    DataT ss[nYdim];
    mGridU2.interpolateUMath(nYdim, knotV, Sv0, Dv0, Sv1, Dv1, v, ss);
    S[0] = ss[0];
  } else {
    mGridU2.interpolateUMath(nYdim, knotV, Sv0, Dv0, Sv1, Dv1, v, S);
  }
  */
  }
};
#endif

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
