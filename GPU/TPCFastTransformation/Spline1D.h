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
/// The spline has n knots x_i. For every knot, the spline value S(x_i) and its derivative S'(x_i) are stored.
/// Inbetween the knots, S(x) is evaluated via interpolation by 3-rd degree polynomials.
///
/// The spline S(x) and its first derivative are continuous at the knots.
/// Depending on the initialization of the derivatives S'(x_i), the second derivative may or may not be continuous.
///
/// --- Knots ---
///
/// The knots are not completely irregular.
/// There is an internal scaled coordinate, called U, where all N knots have integer positions:
/// {U0==0, U1, .., Un-1==Umax}.
/// It is implemented this way for fast matching of any X value to its neighboring knots.
///
/// For example, three knots with U coordinates u_i={0, 3, 5},
/// being stretched on the X segment [0., 1.], will have X coordinates x_i={0., 3./5., 1.}
///
/// For a few reasons, it is better to minimize U-gaps between knots.
/// A spline with knots u_i={0,4,8} is mathematically the same as the spline with knots at u_i={0,1,2},
/// but the later one uses less memory.
///
/// The minimal number of knots is 2.
///
/// --- The output dimensionality ---
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
/// and provide them at each call of the interpolation.
/// To do so, create a spline with nYdimensions=0, create spline parameters for F via SplineHelper1D class,
/// and use special interpolateUMath(..) methods for interpolation.
///
/// This feature allows one to use the same spline object for the approximation of different functions
/// on the same knots.
///
/// ---- Creation of a spline ----
///
/// The splines are best-fit splines. It means, that the spline values S_i and the derivatives D_i at the knots
/// are calibrated such that they minimize the integral difference between S(x) and F(x).
/// This difference is evaluated at all integer values of U coordinate (in particular, at all knots)
/// and at extra nAuxiliaryPoints points between the integers.
///
/// nAuxiliaryPoints can be set as a parameter of approximateFunction() method.
/// With nAuxiliaryPoints==3 the approximation accuracy is noticeably better than the one with 1 or 2.
/// Higher values usually give a little improvement over 3.
///
/// The number of auxiliary points has no influence on the interpolation speed,
/// it can only slow down the approximateFunction() method.
///
/// It is also possible to construct the spline in a classical way - by taking F values only at knots and making
/// the first and the second derivatives of S continuous. To do so, use the corresponding method from SplineHelper1D.
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

///
/// Declare the spline class as a template with one optional parameter
///
template <class DataT, int nYdimT = -1>
class Spline1D;

///
/// The main class specification where the spline dimensionality is set during runtime and read via a class member mYdim
/// DataT is a data type - either double or float.
///
template <class DataT>
class Spline1D<DataT, -1> : public FlatObject
{
 public:
  ///
  /// \brief The struct Knot represents the i-th knot and the segment [knot_i, knot_i+1]
  ///
  struct Knot {
    DataT u;  ///< u coordinate of the knot i (an integer number in float format)
    DataT Li; ///< inverse length of the [knot_i, knot_{i+1}] segment ( == 1./ a (small) integer number)
  };

  /// _____________  Version control __________________________

  /// Version control
  GPUd() static constexpr int getVersion() { return 1; }

    /// _____________  C++ constructors / destructors __________________________

#if !defined(GPUCA_GPUCODE)
  /// Default constructor
  Spline1D() : Spline1D(0, 2) {}

  /// Constructor for a regular spline
  Spline1D(int nYdim, int numberOfKnots);

  /// Constructor for an irregular spline
  Spline1D(int nYdim, int numberOfKnots, const int knotU[]);

  /// Copy constructor
  Spline1D(const Spline1D&);

  /// Assignment operator
  Spline1D& operator=(const Spline1D&);
#else
  /// Disable constructors for the GPU implementation
  Spline1D() CON_DELETE;
  Spline1D(const Spline1D&) CON_DELETE;
  Spline1D& operator=(const Spline1D&) CON_DELETE;
#endif

  /// Destructor
  ~Spline1D() CON_DEFAULT;

  /// _______________  Construction interface  ________________________

#if !defined(GPUCA_GPUCODE)

  /// Constructor for a regular spline.
  void recreate(int nYdim, int numberOfKnots);

  /// Constructor for an irregular spline
  void recreate(int nYdim, int numberOfKnots, const int knotU[]);

  /// approximate a function F with this spline.
  void approximateFunction(double xMin, double xMax,
                           std::function<void(double x, double f[/*mYdim*/])> F,
                           int nAuxiliaryDataPoints = 4);
#endif

  /// _______________  Main functionality   ________________________

  /// Get interpolated value S(x)
  GPUd() void interpolate(DataT x, GPUgeneric() DataT S[/*mYdim*/]) const
  {
    interpolateMath(mYdim, x, S);
  }

  /// Get interpolated value for the first dimension of S(x). (Simplified interface for 1D)
  GPUd() DataT interpolate(DataT x) const
  {
    return interpolateMath(mYdim, x);
  }

    /// _______________  IO   ________________________

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE)
  /// write a class object to the file
  int writeToFile(TFile& outf, const char* name);

  /// read a class object from the file
  static Spline1D* readFromFile(TFile& inpf, const char* name);
#endif

  /// _______________  Getters   ________________________

  /// Get U coordinate of the last knot
  GPUd() int getUmax() const { return mUmax; }

  /// Get number of F dimensions
  GPUd() int getYdimensions() const { return mYdim; }

  /// Get minimal required alignment for the spline parameters
  GPUd() size_t getParameterAlignmentBytes() const
  {
    size_t s = 2 * sizeof(DataT) * mYdim;
    return (s < 16) ? s : 16;
  }

  /// Number of parameters
  GPUd() int getNumberOfParameters() const { return (2 * mYdim) * getNumberOfKnots(); }

  /// Size of the parameter array in bytes
  GPUd() size_t getSizeOfParameters() const { return sizeof(DataT) * getNumberOfParameters(); }

  /// Get a number of knots
  GPUd() int getNumberOfKnots() const { return mNumberOfKnots; }

  /// Get the array of knots
  GPUd() const Knot* getKnots() const { return reinterpret_cast<const Knot*>(mFlatBufferPtr); }

  /// Get i-th knot
  GPUd() const Knot& getKnot(int i) const { return getKnots()[i < 0 ? 0 : (i >= mNumberOfKnots ? mNumberOfKnots - 1 : i)]; }

  /// Get index of an associated knot for a given U coordinate. Performs a boundary check.
  GPUd() int getLeftKnotIndexForU(DataT u) const;

  /// Get u of i-th knot
  GPUd() int getKnotU(int iKnot) const { return (int)(getKnot(iKnot).u + 0.1); }

  /// Get spline parameters
  GPUd() DataT* getParameters() { return mParameters; }

  /// Get spline parameters const
  GPUd() const DataT* getParameters() const { return mParameters; }

  /// _______________  Getters with no boundary check   ________________________

  /// Get i-th knot. No boundary check performed!
  GPUd() const Knot& getKnotNonSafe(int i) const { return getKnots()[i]; }

  /// Get index of an associated knot for a given U coordinate. No bboundaryorder check preformed!
  GPUd() int getLeftKnotIndexForUnonSafe(DataT u) const;

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

  /// Get interpolated value for an mYdim-dimensional S(u) using spline parameters Parameters.
  GPUd() void interpolateU(GPUgeneric() const DataT Parameters[],
                           DataT u, GPUgeneric() DataT S[/*mYdim*/]) const
  {
    interpolateUMath(mYdim, Parameters, u, S);
  }

  /// Same as interpolateU(..) but without boundary checks
  GPUd() void interpolateUnonSafe(GPUgeneric() const DataT Parameters[],
                                  DataT u, GPUgeneric() DataT S[/*mYdim*/]) const
  {
    interpolateUnonSafeMath(mYdim, Parameters, u, S);
  }

  /// The main mathematical utility.
  /// Get interpolated value {S(u): 1D -> mYdim} at the segment [knotL, next knotR]
  /// using the spline values Sl, Sr and the slopes Dl, Dr
  template <typename T>
  GPUd() void interpolateU(const typename Spline1D<DataT>::Knot& knotL,
                           GPUgeneric() const T Sl[/*mYdim*/], GPUgeneric() const T Dl[/*mYdim*/],
                           GPUgeneric() const T Sr[/*mYdim*/], GPUgeneric() const T Dr[/*mYdim*/],
                           DataT u, GPUgeneric() T S[/*mYdim*/]) const
  {
    interpolateUMath(mYdim, knotL, Sl, Dl, Sr, Dr, u, S);
  }

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE) // code invisible on GPU and in the standalone compilation
  /// Test the class functionality
  static int test(const bool draw = 0, const bool drawDataPoints = 1);
#endif

  /// _____________  FlatObject functionality, see FlatObject class for description  ____________

  using FlatObject::getBufferAlignmentBytes;
  using FlatObject::getClassAlignmentBytes;

#if !defined(GPUCA_GPUCODE)
  void cloneFromObject(const Spline1D& obj, char* newFlatBufferPtr);
  void moveBufferTo(char* newBufferPtr);
#endif

  using FlatObject::releaseInternalBuffer;

  void destroy();
  void setActualBufferAddress(char* actualFlatBufferPtr);
  void setFutureBufferAddress(char* futureFlatBufferPtr);

  /// _____________  Interpolation math with external Parameters ____________

  /// Number of parameters
  GPUd() int getNumberOfParametersMath(int nYdim) const { return (2 * nYdim) * getNumberOfKnots(); }

  /// Get interpolated value for an nYdim-dimensional S(u) using spline parameters Parameters.
  GPUd() void interpolateUMath(int nYdim, GPUgeneric() const DataT Parameters[],
                               DataT u, GPUgeneric() DataT S[/*nYdim*/]) const;

  /// Same as interpolateU(..) but with no boundary checks
  GPUd() void interpolateUnonSafeMath(int nYdim, GPUgeneric() const DataT Parameters[],
                                      DataT u, GPUgeneric() DataT S[/*nYdim*/]) const;

  /// The main mathematical utility.
  /// Get interpolated value {S(u): 1D -> mYdim} at the segment [knotL, next knotR]
  /// using the spline values Sl, Sr and the slopes Dl, Dr
  template <typename T>
  GPUd() void interpolateUMath(int nYdim, const typename Spline1D<DataT>::Knot& knotL,
                               GPUgeneric() const T Sl[/*nYdim*/], GPUgeneric() const T Dl[/*nYdim*/],
                               GPUgeneric() const T Sr[/*nYdim*/], GPUgeneric() const T Dr[/*nYdim*/],
                               DataT u, GPUgeneric() T Su[/*nYdim*/]) const;

 protected:
  /// Non-const accessor to knots array
  Knot* getKnots() { return reinterpret_cast<Knot*>(mFlatBufferPtr); }

  /// Non-const accessor to U->knots map
  int* getUtoKnotMap() { return mUtoKnotMap; }

  /// _____________  Interpolation math  ____________

  /// Get interpolated value S(x)
  GPUd() void interpolateMath(int nYdim, DataT x, GPUgeneric() DataT S[/*nYdim*/]) const;

  /// Get interpolated value for the first dimension of S(x). (Simplified interface for 1D)
  GPUd() DataT interpolateMath(int nYdim, DataT x) const;

  /// _____________  Data members  ____________

  int mYdim;          ///< dimentionality of F
  int mNumberOfKnots; ///< n knots on the grid
  int mUmax;          ///< U of the last knot
  DataT mXmin;        ///< X of the first knot
  DataT mXtoUscale;   ///< a scaling factor to convert X to U
  int* mUtoKnotMap;   //! (transient!!) pointer to (integer U -> knot index) map inside the mFlatBufferPtr array
  DataT* mParameters; //! (transient!!) pointer to F-dependent parameters inside the mFlatBufferPtr array
#ifndef GPUCA_ALIROOT_LIB
  ClassDefNV(Spline1D, 1);
#endif
};

///
/// An additional class specification where the spline dimensionality is set at the compile time
/// via a template argument nYdimT
/// It inherits from the main specification above and only replaces some methods by their faster versions.
///
template <class DataT, int nYdimT>
class Spline1D : public Spline1D<DataT, -1>
{
 public:
  typedef Spline1D<DataT, -1> TBase;

  /// _____________  Constructors / destructors __________________________

#if !defined(GPUCA_GPUCODE)
  /// Default constructor
  Spline1D() : Spline1D(2) {}

  /// Constructor for a regular spline && default constructor
  Spline1D(int numberOfKnots) : TBase(nYdimT, numberOfKnots) {}

  /// Constructor for an irregular spline
  Spline1D(int numberOfKnots, const int knotU[]) : TBase(nYdimT, numberOfKnots, knotU) {}

  /// Copy constructor
  Spline1D(const Spline1D& v) : TBase(v) {}

  /// Assignment operator
  Spline1D& operator=(const Spline1D& v) { TBase::operator=(v); };
#else
  /// Disable constructors for the GPU implementation
  Spline1D() CON_DELETE;
  Spline1D(const Spline1D&) CON_DELETE;
  Spline1D& operator=(const Spline1D&) CON_DELETE;
#endif

  /// Destructor
  ~Spline1D() CON_DEFAULT;

#if !defined(GPUCA_GPUCODE)
  /// Constructor for a regular spline.
  void recreate(int numberOfKnots) { TBase::recreate(nYdimT, numberOfKnots); }

  /// Constructor for an irregular spline
  void recreate(int numberOfKnots, const int knotU[]) { TBase::recreate(nYdimT, numberOfKnots, knotU); }
#endif

  /// _______________  Main functionality   ________________________

  /// Get interpolated value S(x)
  GPUd() void interpolate(DataT x, GPUgeneric() DataT S[/*mYdim*/]) const
  {
    interpolateMath(nYdimT, x, S);
  }

  /// Get interpolated value for the first dimension of S(x). (Simplified interface for 1D)
  GPUd() DataT interpolate(DataT x) const
  {
    return interpolateMath(nYdimT, x);
  }

  /// _______________  Getters   ________________________

  /// Get number of F dimensions as constexpr
  GPUd() static constexpr int getYdimensions() { return nYdimT; }

  /// Get minimal required alignment for the spline parameters
  GPUd() static constexpr size_t getParameterAlignmentBytes()
  {
    size_t s = 2 * sizeof(DataT) * nYdimT;
    return (s < 16) ? s : 16;
  }

  /// Number of parameters
  GPUd() int getNumberOfParameters() const { return (2 * nYdimT) * getNumberOfKnots(); }

  /// Size of the parameter array in bytes
  GPUd() size_t getSizeOfParameters() const { return sizeof(DataT) * getNumberOfParameters(); }

  /// ================ Expert tools   ================================

  /// Get interpolated value for an mYdim-dimensional F(u) using spline parameters Parameters.
  GPUd() void interpolateU(GPUgeneric() const DataT Parameters[],
                           DataT u, GPUgeneric() DataT S[/*mYdim*/]) const
  {
    interpolateUMath(nYdimT, Parameters, u, S);
  }

  /// Same as interpolateU(..) but without boundary checks
  GPUd() void interpolateUnonSafe(GPUgeneric() const DataT Parameters[],
                                  DataT u, GPUgeneric() DataT S[/*mYdim*/]) const
  {
    interpolateUnonSafeMath(nYdimT, Parameters, u, S);
  }

  /// The main mathematical utility.
  /// Get interpolated value {S(u): 1D -> mYdim} at the segment [knotL, next knotR]
  /// using the spline values Sl, Sr and the slopes Dl, Dr
  template <typename T>
  GPUd() void interpolateU(const typename Spline1D<DataT>::Knot& knotL,
                           GPUgeneric() const T Sl[/*mYdim*/], GPUgeneric() const T Dl[/*mYdim*/],
                           GPUgeneric() const T Sr[/*mYdim*/], GPUgeneric() const T Dr[/*mYdim*/],
                           DataT u, GPUgeneric() T S[/*mYdim*/]) const
  {
    interpolateUMath(nYdimT, knotL, Sl, Dl, Sr, Dr, u, S);
  }

    /// _______________  IO   ________________________

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE)
  /// read a class object from the file
  static Spline1D* readFromFile(TFile& inpf, const char* name)
  {
    return reinterpret_cast<Spline1D*>(TBase::readFromFile(inpf, name));
  }
#endif

  /// _______________  Suppress some base class methods   ________________________
 private:
#if !defined(GPUCA_GPUCODE)
  void recreate(int nYdim, int numberOfKnots){};
  void recreate(int nYdim, int numberOfKnots, const int knotU[]){};
#endif
 public:
  using TBase::getNumberOfKnots;
  using TBase::interpolateMath;
  using TBase::interpolateUMath;
  using TBase::interpolateUnonSafeMath;
#ifndef GPUCA_ALIROOT_LIB
  //ClassDefNV(Spline1D, 1);
#endif
};

///
/// ========================================================================================================
///       Inline implementations
/// ========================================================================================================
///

#if !defined(GPUCA_GPUCODE)

template <class DataT>
GPUdi() Spline1D<DataT>::Spline1D(int nYdim, int numberOfKnots)
  : FlatObject()
{
  recreate(nYdim, numberOfKnots);
}

template <class DataT>
GPUdi() Spline1D<DataT>::Spline1D(int nYdim, int numberOfKnots, const int knotU[])
  : FlatObject()
{
  recreate(nYdim, numberOfKnots, knotU);
}

template <class DataT>
GPUdi() Spline1D<DataT>::Spline1D(const Spline1D& spline)
  : FlatObject()
{
  cloneFromObject(spline, nullptr);
}

template <class DataT>
GPUd() Spline1D<DataT>& Spline1D<DataT>::operator=(const Spline1D<DataT>& spline)
{
  cloneFromObject(spline, nullptr);
  return *this;
}
#endif

template <class DataT>
GPUdi() int Spline1D<DataT>::getLeftKnotIndexForU(DataT u) const
{
  /// Get i: u is in [knot_i, knot_{i+1}) segments
  /// when u is otside of [0, mUmax], return a corresponding edge segment
  int iu = (int)u;
  if (iu < 0) {
    iu = 0;
  }
  if (iu > mUmax) {
    iu = mUmax;
  }
  return getUtoKnotMap()[iu];
}

template <class DataT>
GPUdi() int Spline1D<DataT>::getLeftKnotIndexForUnonSafe(DataT u) const
{
  /// Get i: u is in [knot_i, knot_{i+1}) segment
  /// no boundary check! u must be in [0,mUmax]
  return getUtoKnotMap()[(int)u];
}

template <class DataT>
GPUdi() void Spline1D<DataT>::setXrange(DataT xMin, DataT xMax)
{
  mXmin = xMin;
  double l = ((double)xMax) - xMin;
  if (l < 1.e-8) {
    l = 1.e-8;
  }
  mXtoUscale = mUmax / l;
}

template <class DataT>
GPUdi() void Spline1D<DataT>::interpolateMath(int nYdim, DataT x, GPUgeneric() DataT S[]) const
{
  /// Get interpolated value S(x)
  interpolateUMath(nYdim, mParameters, convXtoU(x), S);
}

template <class DataT>
GPUdi() DataT Spline1D<DataT>::interpolateMath(int nYdim, DataT x) const
{
  /// Simplified interface for 1D: get interpolated value for the first dimension of S(x)
  DataT u = convXtoU(x);
  int iknot = getLeftKnotIndexForU(u);
  const DataT* d = mParameters + (2 * nYdim) * iknot;
  DataT S = 0.;
  interpolateUMath(nYdim, getKnotNonSafe(iknot), &(d[0]), &(d[nYdim]),
                   &(d[2 * nYdim]), &(d[3 * nYdim]), u, &S);
  return S;
}

template <class DataT>
GPUdi() void Spline1D<DataT>::interpolateUMath(int nYdim, GPUgeneric() const DataT parameters[], DataT u, GPUgeneric() DataT S[/*mYdim*/]) const
{
  /// Get interpolated value S(u) using given spline parameters with a boundary check
  int iknot = getLeftKnotIndexForU(u);
  const DataT* d = parameters + (2 * nYdim) * iknot;
  interpolateUMath(nYdim, getKnotNonSafe(iknot), &(d[0]), &(d[nYdim]), &(d[2 * nYdim]), &(d[3 * nYdim]), u, S);
}

template <class DataT>
GPUdi() void Spline1D<DataT>::interpolateUnonSafeMath(int nYdim, GPUgeneric() const DataT parameters[], DataT u, GPUgeneric() DataT S[/*mYdim*/]) const
{
  /// Get interpolated value S(u) using given spline parameters without boundary check
  int iknot = getLeftKnotIndexForUnonSafe(u);
  const DataT* d = parameters + (2 * nYdim) * iknot;
  interpolateUMath(nYdim, getKnotNonSafe(iknot), &(d[0]), &(d[nYdim]), &(d[2 * nYdim]), &(d[3 * nYdim]), u, S);
}

template <class DataT>
template <typename T>
GPUdi() void Spline1D<DataT>::interpolateUMath(int nYdim, const typename Spline1D<DataT>::Knot& knotL,
                                               GPUgeneric() const T Sl[/*nYdim*/], GPUgeneric() const T Dl[/*nYdim*/],
                                               GPUgeneric() const T Sr[/*nYdim*/], GPUgeneric() const T Dr[/*nYdim*/],
                                               DataT u, GPUgeneric() T S[/*nYdim*/]) const
{
  /// A static method.
  /// Gives interpolated value of N-dimensional S(u) at u
  /// input: Sl,Dl,Sr,Dr[nYdim] - N-dim function values and slopes at knots {knotL,knotR}
  /// output: S[nYdim] - N-dim interpolated value for S(u)

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

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
