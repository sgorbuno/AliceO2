// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  Spline1D.h
/// \brief Definition of Spline1D class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#ifndef ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SPLINE1D_H
#define ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SPLINE1D_H

#include "GPUCommonDef.h"
#include "FlatObject.h"
#if !defined(GPUCA_GPUCODE)
#include <functional>
#endif

class TFile;

namespace GPUCA_NAMESPACE
{
namespace gpu
{
/// The Spline1D class performs a cubic spline interpolation on a one-dimensional non-uniform grid.
///
/// The class is a flat C structure. It inherits from the FlatObject.
/// No virtual methods, no ROOT types are used.
///
/// --- Interpolation ---
///
/// The spline S(x) approximates a function F(x):[Xmin,Xmax]->Y
/// X is one-dimensional, Y may be multi-dimensional.
///
/// The spline has n knots x_i. The spline value S(x_i) and its derivative S'(x_i) are stored for every knot.
/// Inbetween the knots, S(x) is evaluated via interpolation by 3-rd degree polynomials.
///
/// The spline S(x) and its first derivative are continuous.
/// Depending on the initialization of the derivatives S'(x_i),
/// the second derivative may or may not be continuous at the knots.
///
/// --- Knots ---
///
/// The knots are not entirely irregular.
/// There is an internal scaled coordinate, called U, where all N knots have some integer positions:
/// {U0==0, U1, .., Un-1==Umax}.
/// It is implemented this way for fast matching of any X value to its neighboring knots.
///
/// For example, three knots with U coordinates u_i={0, 3, 5},
/// being stretched on the X segment [0., 1.], will have X coordinates x_i={0., 3./5., 1.}
///
/// For a few reasons, it is better to keep U-gaps between the knots minimal.
/// A spline with knots u_i={0,1,2} is mathematically the same as the spline with knots u_i={0,2,4}.
/// However, it uses less memory.
///
/// The minimal number of knots is 2.
///
/// --- Output dimensionality ---
///
/// There are two ways to set the dimensionality of Y - either in the constructor or as a template argument:
///
/// Spline1D<float> s( nYdimensions, nKnots );
/// Spline1D<float, nYdimensions> s( nKnots );
///
/// The second implementation works faster. Use it when nYdimensions is known at the compile time.
///
/// ---- External storage of spline parameters for a given F ---
///
/// One can store all F-dependent spline parameters outside of the spline object
/// and provide them at each interpolation call.
/// To do so, create a spline with nYdimensions=0; create spline parameters for F via SplineHelper1D class;
/// then use special interpolateUMath(..) methods for interpolation.
///
/// This feature allows one to use the same spline object for the approximation of different functions
/// on the same grid of knots.
///
/// ---- Creation of a spline ----
///
/// The spline is supposed to be a best-fit spline, created by the approximateFunction() method.
///
/// Best-fit means that the spline values S_i and its derivatives D_i at the knots
/// are adjusted to minimize the overall difference between S(x) and F(x).
/// The spline constructed this way is much more accurate than a classical interpolation spline.
///
/// The difference to F() is minimized at all integer values of U coordinate (in particular, at all knots)
/// and at extra nAuxiliaryPoints points between the integer numbers.
///
/// nAuxiliaryPoints is given as a parameter of approximateFunction() method.
/// With nAuxiliaryPoints==3, the approximation accuracy is noticeably better than the one with 1 or 2.
/// Higher values usually give a little improvement over 3.
///
/// The number of auxiliary points does not influence the interpolation speed,
/// but a high number can slow down the spline's creation.
///
/// It is also possible to construct the spline classically - by taking F(x) values only at knots and making
/// the first and the second derivatives of S(x) continuous. Use the corresponding method from SplineHelper1D.
///
/// ---- Example of creating a spline ----
///
///  auto F = [&](double x, double &f) { // a function to be approximated
///   f[0] = x*x+3.f; // F(x)
///  };
///
///  const int nKnots = 3;
///
///  int knots[nKnots] = {0, 1, 5}; // relative(!) knot positions
///
///  Spline1D<float,1> spline( nKnots, knots ); // create 1-dimensional spline with the knots
///
///  spline.approximateFunction(0., 1., F); // let the spline approximate F on a segment [0., 1.]
///
///  float s = spline.interpolate(0.2); // interpolated value at x==0.2
///
///  --- See also Spline1D::test() method for examples
///
///

// ==================================================================================================

/// The class Spline1DContainer is a base class of Spline1D.
/// It contains all the class members and methods which only depends on the DataT data type
/// and do not depend on other template paramaters of Spline1D.
/// It also contaions all non-inined methods with the implementation in .cxx file.
///
/// DataT is a data type, which is supposed to be either double or float.
/// For other possible data types one has to add the corresponding instantiation line
/// at the end of the .cxx file
///
template <typename DataT>
class Spline1DContainer : public FlatObject
{
 public:
  /// A named enumeration for the safety level used by some methods
  enum SafetyLevel { kNotSafe,
                     kSafe };

  /// The struct Knot represents the i-th knot and the segment [knot_i, knot_i+1]
  ///
  struct Knot {
    DataT u;  ///< u coordinate of the knot i (an integer number in float format)
    DataT Li; ///< inverse length of the [knot_i, knot_{i+1}] segment ( == 1./ a (small) integer )
    /// Get u as an integer
    GPUd() int getU() const { return (int)(u + 0.1); }
  };

  /// _____________  Version control __________________________

  /// Version control
  GPUd() static constexpr int getVersion() { return 1; }

    /// _____________  C++ constructors / destructors __________________________

#if !defined(GPUCA_GPUCODE)
  /// Default constructor, needed by the Root IO
  Spline1DContainer() : FlatObject() { recreate(0, 2); }

 protected:
  /// Copy constructor
  Spline1DContainer(const Spline1DContainer&);
  /// Assignment operator
  Spline1DContainer& operator=(const Spline1DContainer& v);
#else
  /// Disable constructors for the GPU implementation
  Spline1DContainer() CON_DELETE;
  Spline1DContainer(const Spline1DContainer&) CON_DELETE;
#endif

 public:
  /// Destructor
  ~Spline1DContainer() CON_DEFAULT;

  /// _______________  Construction interface  ________________________

#if !defined(GPUCA_GPUCODE)
  /// Constructor for a regular spline
  void recreate(int nYdim, int numberOfKnots);

  /// Constructor for an irregular spline
  void recreate(int nYdim, int numberOfKnots, const int knotU[]);

  /// approximate a function F with this spline
  void approximateFunction(double xMin, double xMax,
                           std::function<void(double x, double f[/*mYdim*/])> F,
                           int nAuxiliaryDataPoints = 4);
#endif

  /// _______________  IO   ________________________

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE)
  /// write a class object to the file
  int writeToFile(TFile& outf, const char* name);

  /// read a class object from the file
  static Spline1DContainer* readFromFile(TFile& inpf, const char* name);
#endif

  /// _______________  Getters   ________________________

  /// Get U coordinate of the last knot
  GPUd() int getUmax() const { return mUmax; }

  /// Get number of Y dimensions
  GPUd() int getYdimensions() const { return mYdim; }

  /// Get minimal required alignment for the spline parameters
  GPUd() size_t getParameterAlignmentBytes() const
  {
    size_t s = 2 * sizeof(DataT) * mYdim;
    return (s < 16) ? s : 16;
  }

  /// Number of parameters
  GPUd() int getNumberOfParameters() const { return calcNumberOfParameters(mYdim); }

  /// Size of the parameter array in bytes
  GPUd() size_t getSizeOfParameters() const { return sizeof(DataT) * getNumberOfParameters(); }

  /// Get a number of knots
  GPUd() int getNumberOfKnots() const { return mNumberOfKnots; }

  /// Get the array of knots
  GPUd() const Knot* getKnots() const { return reinterpret_cast<const Knot*>(mFlatBufferPtr); }

  /// Get i-th knot
  template <SafetyLevel SafeT = SafetyLevel::kSafe>
  GPUd() const Knot& getKnot(int i) const
  {
    if (SafeT == SafetyLevel::kSafe) {
      i = (i < 0) ? 0 : (i >= mNumberOfKnots ? mNumberOfKnots - 1 : i);
    }
    return getKnots()[i];
  }

  /// Get index of an associated knot for a given U coordinate. Performs a boundary check.
  template <SafetyLevel SafeT = SafetyLevel::kSafe>
  GPUd() int getLeftKnotIndexForU(DataT u) const;

  /// Get spline parameters
  GPUd() DataT* getParameters() { return mParameters; }

  /// Get spline parameters const
  GPUd() const DataT* getParameters() const { return mParameters; }

  /// _______________  Technical stuff  ________________________

  /// Get a map (integer U -> corresponding knot index)
  GPUd() const int* getUtoKnotMap() const { return mUtoKnotMap; }

  /// Convert X coordinate to U
  GPUd() DataT convXtoU(DataT x) const { return (x - mXmin) * mXtoUscale; }

  /// Convert U coordinate to X
  GPUd() DataT convUtoX(DataT u) const { return mXmin + u / mXtoUscale; }

  /// Get Xmin
  GPUd() DataT getXmin() const { return mXmin; }

  /// Get XtoUscale
  GPUd() DataT getXtoUscale() const { return mXtoUscale; }

  /// Set X range
  void setXrange(DataT xMin, DataT xMax);

  /// Print method
  void print() const;

  ///  _______________  Expert tools  _______________

  /// Number of parameters for a given Y dimensions
  GPUd() int calcNumberOfParameters(int nYdim) const { return (2 * nYdim) * getNumberOfKnots(); }

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE) // code invisible on GPU and in the standalone compilation
  /// Test the class functionality
  static int test(const bool draw = 0, const bool drawDataPoints = 1);
#endif

  /// _____________  FlatObject functionality, see FlatObject class for description  ____________

  using FlatObject::getBufferAlignmentBytes;
  using FlatObject::getClassAlignmentBytes;

#if !defined(GPUCA_GPUCODE)
  void cloneFromObject(const Spline1DContainer& obj, char* newFlatBufferPtr);
  void moveBufferTo(char* newBufferPtr);
#endif

  using FlatObject::releaseInternalBuffer;

  void destroy();
  void setActualBufferAddress(char* actualFlatBufferPtr);
  void setFutureBufferAddress(char* futureFlatBufferPtr);

 protected:
  /// Non-const accessor to the knots array
  Knot* getKnots() { return reinterpret_cast<Knot*>(mFlatBufferPtr); }

  /// Non-const accessor to U->knots map
  int* getUtoKnotMap() { return mUtoKnotMap; }

  /// _____________  Data members  ____________

  int mYdim;          ///< dimentionality of F
  int mNumberOfKnots; ///< n knots on the grid
  int mUmax;          ///< U of the last knot
  DataT mXmin;        ///< X of the first knot
  DataT mXtoUscale;   ///< a scaling factor to convert X to U
  int* mUtoKnotMap;   //! (transient!!) pointer to (integer U -> knot index) map inside the mFlatBufferPtr array
  DataT* mParameters; //! (transient!!) pointer to F-dependent parameters inside the mFlatBufferPtr array
#ifndef GPUCA_ALIROOT_LIB
  ClassDefNV(Spline1DContainer, 1);
#endif
};

#if !defined(GPUCA_GPUCODE)
template <typename DataT>
GPUdi() Spline1DContainer<DataT>::Spline1DContainer(const Spline1DContainer<DataT>& spline)
  : FlatObject()
{
  cloneFromObject(spline, nullptr);
}
template <typename DataT>
GPUdi() Spline1DContainer<DataT>& Spline1DContainer<DataT>::operator=(const Spline1DContainer<DataT>& spline)
{
  cloneFromObject(spline, nullptr);
  return *this;
}
#endif

template <typename DataT>
template <typename Spline1DContainer<DataT>::SafetyLevel SafeT>
GPUdi() int Spline1DContainer<DataT>::getLeftKnotIndexForU(DataT u) const
{
  /// Get i: u is in [knot_i, knot_{i+1}) segments
  /// when u is otside of [0, mUmax], return a corresponding edge segment
  int iu = (int)u;
  if (SafeT == SafetyLevel::kSafe) {
    iu = (iu < 0) ? 0 : (iu > mUmax ? mUmax : iu);
  }
  return getUtoKnotMap()[iu];
}

template <typename DataT>
GPUdi() void Spline1DContainer<DataT>::setXrange(DataT xMin, DataT xMax)
{
  mXmin = xMin;
  double l = ((double)xMax) - xMin;
  if (l < 1.e-8) {
    l = 1.e-8;
  }
  mXtoUscale = mUmax / l;
}

/// ==================================================================================================
///
/// The class declares methods which are common for all Spline1D specifications.
/// Implementation of the methods may differ depending on the template arguments.
///
template <class DataT, int nYdimT, bool useYdimT, bool fixedMemAllocT>
class Spline1DCommon
  : public Spline1DContainer<DataT>
{
  typedef Spline1DContainer<DataT> TBase;
  typedef typename TBase::SafetyLevel SafetyLevel;
  typedef typename TBase::Knot Knot;

 public:
  /// _____________  C++ constructors / destructors __________________________

#if !defined(GPUCA_GPUCODE)
 protected:
  /// Default constructor
  Spline1DCommon() : TBase() {}
  /// Copy constructor
  Spline1DCommon(const Spline1DCommon& v) : TBase(v){};
#else
  /// Disable constructors for the GPU implementation
  Spline1DCommon() CON_DELETE;
  Spline1DCommon(const Spline1DCommon&) CON_DELETE;
#endif

 public:
  /// Destructor
  ~Spline1DCommon() CON_DEFAULT;

  /// _______________  Main functionality   ________________________

  /// Get interpolated value S(x)
  GPUd() void interpolate(DataT x, GPUgeneric() DataT S[/*mYdim*/]) const
  {
    interpolateU<SafetyLevel::kSafe>(mYdim, mParameters, convXtoU(x), S);
  }

  /// Get interpolated value for the first dimension of S(x). (Simplified interface for 1D)
  GPUd() DataT interpolate(DataT x) const;

  using TBase::convXtoU;
  using TBase::getKnot;
  using TBase::getKnots;
  using TBase::getNumberOfKnots;

 protected:
  using TBase::mParameters;
  using TBase::mYdim;

  /// Get interpolated value for an nYdim-dimensional S(u) using spline parameters Parameters.
  template <SafetyLevel SafeT = SafetyLevel::kSafe>
  GPUd() void interpolateU(int nYdim, GPUgeneric() const DataT Parameters[],
                           DataT u, GPUgeneric() DataT S[/*nYdim*/]) const;

  /// The main mathematical utility.
  /// Get interpolated value {S(u): 1D -> nYdim} at the segment [knotL, next knotR]
  /// using the spline values Sl, Sr and the slopes Dl, Dr
  template <typename T>
  GPUd() void interpolateU(int nYdim, const Knot& knotL,
                           GPUgeneric() const T Sl[/*mYdim*/], GPUgeneric() const T Dl[/*mYdim*/],
                           GPUgeneric() const T Sr[/*mYdim*/], GPUgeneric() const T Dr[/*mYdim*/],
                           DataT u, GPUgeneric() T S[/*mYdim*/]) const;

#ifndef GPUCA_ALIROOT_LIB
  ClassDefNV(Spline1DCommon, 1);
#endif
};

template <typename DataT, int nYdimT, bool useYdimT, bool fixedMemAllocT>
GPUdi() DataT Spline1DCommon<DataT, nYdimT, useYdimT, fixedMemAllocT>::interpolate(DataT x) const
{
  /// Simplified interface for 1D: get interpolated value for the first dimension of S(x)
  int nYdim = mYdim;
  DataT u = convXtoU(x);
  int iknot = TBase::template getLeftKnotIndexForU<SafetyLevel::kSafe>(u);
  const DataT* d = mParameters + (2 * nYdim) * iknot;
  DataT S = 0.;
  ((Spline1DCommon<DataT, 1, true, true>*)this)->interpolateU(nYdim, getKnots()[iknot], &(d[0]), &(d[nYdim]), &(d[2 * nYdim]), &(d[3 * nYdim]), u, &S);
  return S;
}

template <typename DataT, int nYdimT, bool useYdimT, bool fixedMemAllocT>
template <typename Spline1DContainer<DataT>::SafetyLevel SafeT>
GPUdi() void Spline1DCommon<DataT, nYdimT, useYdimT, fixedMemAllocT>::interpolateU(
  int nYdim, GPUgeneric() const DataT parameters[], DataT u, GPUgeneric() DataT S[/*mYdim*/]) const
{
  /// Get interpolated value S(u) using given spline parameters with a boundary check
  //int nYdim = mYdim;
  int iknot = TBase::template getLeftKnotIndexForU<SafeT>(u);
  const DataT* d = parameters + (2 * nYdim) * iknot;
  interpolateU(nYdim, getKnots()[iknot], &(d[0]), &(d[nYdim]), &(d[2 * nYdim]), &(d[3 * nYdim]), u, S);
}

template <typename DataT, int nYdimT, bool useYdimT, bool fixedMemAllocT>
template <typename T>
GPUdi() void Spline1DCommon<DataT, nYdimT, useYdimT, fixedMemAllocT>::interpolateU(
  int nYdim, const typename Spline1DContainer<DataT>::Knot& knotL,
  GPUgeneric() const T Sl[/*nYdim*/], GPUgeneric() const T Dl[/*nYdim*/],
  GPUgeneric() const T Sr[/*nYdim*/], GPUgeneric() const T Dr[/*nYdim*/],
  DataT u, GPUgeneric() T S[/*nYdim*/]) const
{
  /// A static method.
  /// Gives interpolated value of N-dimensional S(u) at u
  /// input: Sl,Dl,Sr,Dr[nYdim] - N-dim function values and slopes at knots {knotL,knotR}
  /// output: S[nYdim] - N-dim interpolated value for S(u)

  //int nYdim = mYdim;
  T uu = T(u - knotL.u);
  T li = T(knotL.Li);
  T v = uu * li; // scaled u
  for (int dim = 0; dim < nYdim; ++dim) {
    T df = (Sr[dim] - Sl[dim]) * li;
    T a = Dl[dim] + Dr[dim] - df - df;
    T b = df - Dl[dim] - a;
    S[dim] = ((a * v + b) * v + Dl[dim]) * uu + Sl[dim];
  }

  /* another way to calculate f(u):
  T uu = T(u - knotL.u);
  T v = uu * T(knotL.Li); // scaled u
  T vm1 = v-1;
  T v2 = v * v;
  float cSr = v2*(3-2*v);
  float cSl = 1-cSr;
  float cDl = v*vm1*vm1*knotL.L;
  float cDr = v2*vm1*knotL.L;
  return cSl*Sl + cSr*Sr + cDl*Dl + cDr*Dr;
  */
}

/// ==================================================================================================
///
/// The class declares methods which are different for different Spline1D specifications
///
template <typename DataT, int nYdimT, bool useYdimT, bool fixedMemAllocT>
class Spline1DSpecific;

/// ==================================================================================================
///
/// Specification where nYdim is taken during runtime as a constructor argument
///
template <typename DataT, int nYdimT, bool fixedMemAllocT>
class Spline1DSpecific<DataT, nYdimT, false, fixedMemAllocT>
  : public Spline1DCommon<DataT, nYdimT, false, fixedMemAllocT>
{
  typedef Spline1DContainer<DataT> TVeryBase;
  typedef Spline1DCommon<DataT, nYdimT, false, fixedMemAllocT> TBase;

 public:
  typedef typename TVeryBase::SafetyLevel SafetyLevel;
  typedef typename TVeryBase::Knot Knot;

  /// _____________  C++ constructors / destructors __________________________

#if !defined(GPUCA_GPUCODE)
  /// Default constructor
  Spline1DSpecific() : Spline1DSpecific(0, 2) {}

  /// Constructor for a regular spline
  Spline1DSpecific(int nYdim, int numberOfKnots) : TBase()
  {
    TBase::recreate(nYdim, numberOfKnots);
  }
  /// Constructor for an irregular spline
  Spline1DSpecific(int nYdim, int numberOfKnots, const int knotU[]) : TBase()
  {
    TBase::recreate(nYdim, numberOfKnots, knotU);
  }
  /// Copy constructor
  Spline1DSpecific(const Spline1DSpecific& v) : TBase(v) {}
#else
  /// Disable constructors for the GPU implementation
  Spline1DSpecific() CON_DELETE;
  Spline1DSpecific(const Spline1DSpecific&) CON_DELETE;
#endif

  /// _____________  C++ constructors / destructors __________________________

  /// Destructor
  ~Spline1DSpecific() CON_DEFAULT;

  ///  _______  Expert tools: interpolation with given nYdim and external Parameters _______

  /// Get interpolated value for an nYdim-dimensional S(u) using spline parameters Parameters.
  template <SafetyLevel SafeT = SafetyLevel::kSafe>
  GPUd() void interpolateU(int nYdim, GPUgeneric() const DataT Parameters[],
                           DataT u, GPUgeneric() DataT S[/*nYdim*/]) const
  {
    TBase::template interpolateU<SafeT>(nYdim, Parameters, u, S);
  }

 public:
#ifndef GPUCA_ALIROOT_LIB
  ClassDefNV(Spline1DSpecific, 1);
#endif
};

/// ==================================================================================================
///
/// Specification where nYdim is taken from a template argument nYdimT at the compile time
///
template <typename DataT, int nYdimT, bool fixedMemAllocT>
class Spline1DSpecific<DataT, nYdimT, true, fixedMemAllocT>
  : public Spline1DCommon<DataT, nYdimT, true, fixedMemAllocT>
{
  typedef Spline1DContainer<DataT> TVeryBase;
  typedef Spline1DCommon<DataT, nYdimT, true, fixedMemAllocT> TBase;

 public:
  typedef typename TVeryBase::SafetyLevel SafetyLevel;
  typedef typename TVeryBase::Knot Knot;

  /// _____________  C++ constructors / destructors __________________________

#if !defined(GPUCA_GPUCODE)
  /// Default constructor
  Spline1DSpecific() : Spline1DSpecific(2) {}

  /// Constructor for a regular spline
  Spline1DSpecific(int numberOfKnots) : TBase()
  {
    recreate(numberOfKnots);
  }
  /// Constructor for an irregular spline
  Spline1DSpecific(int numberOfKnots, const int knotU[])
    : TBase()
  {
    recreate(numberOfKnots, knotU);
  }
  /// Copy constructor
  Spline1DSpecific(const Spline1DSpecific& v) : TBase(v){};
#else
  /// Disable constructors for the GPU implementation
  Spline1DSpecific() CON_DELETE;
  Spline1DSpecific(const Spline1DSpecific&) CON_DELETE;
#endif

  /// Destructor
  ~Spline1DSpecific() CON_DEFAULT;

  /// _______________  Construction interface  ________________________

#if !defined(GPUCA_GPUCODE)

  /// Constructor for a regular spline
  void recreate(int numberOfKnots) { TBase::recreate(nYdimT, numberOfKnots); }

  /// Constructor for an irregular spline
  void recreate(int numberOfKnots, const int knotU[])
  {
    TBase::recreate(nYdimT, numberOfKnots, knotU);
  }
#endif

  /// _______________  Getters   ________________________

  /// Get number of Y dimensions
  GPUd() constexpr int getYdimensions() const { return nYdimT; }

  /// Get minimal required alignment for the spline parameters
  GPUd() constexpr size_t getParameterAlignmentBytes() const
  {
    size_t s = 2 * sizeof(DataT) * nYdimT;
    return (s < 16) ? s : 16;
  }

  /// Number of parameters
  GPUd() int getNumberOfParameters() const { return (2 * nYdimT) * getNumberOfKnots(); }

  /// Size of the parameter array in bytes
  GPUd() size_t getSizeOfParameters() const { return (sizeof(DataT) * 2 * nYdimT) * getNumberOfKnots(); }

  ///  _______  Expert tools: interpolation with given nYdim and external Parameters _______

  /// Get interpolated value for an nYdimT-dimensional S(u) using spline parameters Parameters.
  template <SafetyLevel SafeT = SafetyLevel::kSafe>
  GPUd() void interpolateU(GPUgeneric() const DataT Parameters[],
                           DataT u, GPUgeneric() DataT S[/*nYdim*/]) const
  {
    TBase::template interpolateU<SafeT>(nYdimT, Parameters, u, S);
  }

  using TBase::getNumberOfKnots;

  /// _______________  Suppress some base class methods   ________________________
 private:
#if !defined(GPUCA_GPUCODE)
  using TBase::recreate;
#endif

 public:
#ifndef GPUCA_ALIROOT_LIB
  ClassDefNV(Spline1DSpecific, 1);
#endif
};

/// ==================================================================================================
///
/// Declare the Spline1D class as a template with one optional parameter
///
template <typename DataT, int nYdimT = -666>
class Spline1D;

/// ==================================================================================================
///
/// Spline1D specification with one template argument.
/// nYdim is given at the runtime via constructor,
/// nYdim-dependent temporary memory might be dynamically allocated in the stack during calculations
/// (not relevant for the actual code)
///
/// \param DataT data type: float or double
///
template <typename DataT>
class Spline1D<DataT, -666> : public Spline1DSpecific<DataT, 0, false, false>
{
  typedef Spline1DContainer<DataT> TVeryBase;
  typedef Spline1DSpecific<DataT, 0, false, false> TBase;

 public:
  typedef typename TVeryBase::SafetyLevel SafetyLevel;
  typedef typename TVeryBase::Knot Knot;

  using TBase::TBase; // inherit constructors

#if !defined(GPUCA_GPUCODE)
  /// Assignment operator
  Spline1D& operator=(const Spline1D& v) { return (Spline1D&)TVeryBase::operator=(v); }
#endif

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE)
  /// read a class object from the file
  static Spline1D* readFromFile(TFile& inpf, const char* name)
  {
    return (Spline1D*)TVeryBase::readFromFile(inpf, name);
  }
#endif

  ClassDefNV(Spline1D, 0);
};

/// ==================================================================================================
///
/// Spline1D specification with two template arguments.
///
/// Case 1: nYdimT >= 0
///    nYdim is set at the compile time via nYdimT.
/// Case 2: nYdimT < 0
///    nYdim is given at the runtime via constructor,
///    nYdim-dependent temporary memory for calculations is statically allocated.
///    abs(nYdimT) is the maximal value of nYdim.
///    (not relevant for the actual code)
///
/// \param DataT data type: float or double
/// \param nYdimT >=0: number of Y dimensions, < 0 : max possible number of Y dimensions
///

template <typename DataT, int nYdimT>
class Spline1D : public Spline1DSpecific<DataT, nYdimT, (nYdimT >= 0), true>
{
 public:
  typedef Spline1DContainer<DataT> TVeryBase;
  typedef Spline1DSpecific<DataT, nYdimT, (nYdimT >= 0), true> TBase;

  typedef typename TVeryBase::SafetyLevel SafetyLevel;
  typedef typename TVeryBase::Knot Knot;

  using TBase::TBase; // inherit constructors

#if !defined(GPUCA_GPUCODE)
  /// Assignment operator
  Spline1D& operator=(const Spline1D& v) { return (Spline1D&)TVeryBase::operator=(v); }
#endif

  /// cast operator
  operator Spline1D<DataT>&() { return *((Spline1D<DataT>*)this); }

  /// const cast operator
  operator const Spline1D<DataT>&() const { return *((const Spline1D<DataT>*)this); }

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE)
  /// read a class object from the file
  static Spline1D* readFromFile(TFile& inpf, const char* name)
  {
    return (Spline1D*)TVeryBase::readFromFile(inpf, name);
  }
#endif

  ClassDefNV(Spline1D, 0);
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
