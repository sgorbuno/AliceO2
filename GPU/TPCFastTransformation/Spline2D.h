// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  Spline2D.h
/// \brief Definition of Spline2D class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#ifndef ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SPLINE2D_H
#define ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SPLINE2D_H

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
///
/// The Spline2D class performs a cubic spline interpolation on an two-dimensional nonunifom grid.
/// The class is an extension of the Spline1D class.
/// See Spline1D.h for more details.
///
/// The spline S(x1,x2) approximates a function F(x1,x2):R^2->R^m,
/// with 2-dimensional domain and multi-dimensional codomain.
/// x1,x2 belong to [x1min,x1max] x [x2min,x2max].
///
/// --- Example of creating a spline ---
///
///  auto F = [&](float x1, float x2, float f[] ) {
///   f[0] = 1.f + x1 + x2*x2; // F(x1,x2)
///  };
///  const int nKnotsU=2;
///  const int nKnotsV=3;
///  int knotsU[nKnotsU] = {0, 1};
///  int knotsV[nKnotsV] = {0, 2, 5};
///  Spline2D<float,1> spline(nKnotsU, knotsU, nKnotsV, knotsV ); // spline with 1-dimensional codomain
///  spline.approximateFunction(0., 1., 0.,1., F); //initialize spline to approximate F on [0., 1.]x[0., 1.] area
///  float S = spline.interpolate(.1, .3 ); // interpolated value at (.1,.3)
///
///  --- See also Spline2D::test();
///

/// Base class to store data members and non-inline methods
template <typename DataT>
class Spline2DBase : public FlatObject
{
 public:
  typedef typename Spline1D<DataT>::SafeFlag SafeFlag;
  typedef typename Spline1D<DataT>::Knot Knot;

  /// _____________  Version control __________________________

  /// Version control
  GPUd() static constexpr int getVersion() { return (1 << 16) + Spline1D<DataT>::getVersion(); }

  /// _____________  C++ constructors / destructors __________________________

#if !defined(GPUCA_GPUCODE)
  /// Default constructor
  Spline2DBase() : Spline2DBase(0, 2, 2) {}

  /// Constructor for a regular spline
  Spline2DBase(int nYdim, int numberOfKnotsU1, int numberOfKnotsU2);

  /// Constructor for an irregular spline
  Spline2DBase(int nYdim, int numberOfKnotsU1, const int knotU1[], int numberOfKnotsU2, const int knotU2[]);

  /// Copy constructor
  Spline2DBase(const Spline2DBase&);

  /// Assignment operator
  Spline2DBase& operator=(const Spline2DBase&);
#else
  /// Disable constructors for the GPU implementation
  Spline2DBase() CON_DELETE;
  Spline2DBase(const Spline2DBase&) CON_DELETE;
  Spline2DBase& operator=(const Spline2DBase&) CON_DELETE;
#endif

  /// Destructor
  ~Spline2DBase() CON_DEFAULT;

  /// _______________  Construction interface  ________________________

#if !defined(GPUCA_GPUCODE)

  /// Constructor for a regular spline
  void recreate(int nYdim, int numberOfKnotsU1, int numberOfKnotsU2);

  /// Constructor for an irregular spline
  void recreate(int nYdim, int numberOfKnotsU1, const int knotU1[], int numberOfKnotsU2, const int knotU2[]);

  /// approximate a function F with this spline
  void approximateFunction(double x1Min, double x1Max, double x2Min, double x2Max,
                           std::function<void(double x1, double x2, double f[/*mYdim*/])> F,
                           int nAuxiliaryDataPointsU1 = 4, int nAuxiliaryDataPointsU2 = 4);
#endif

  /// _______________  Main functionality   ________________________

  /// Get interpolated value S(x)
  GPUd() void interpolate(DataT x1, DataT x2, GPUgeneric() DataT S[/*mYdim*/]) const
  {
    interpolateMath(mYdim, x1, x2, S);
  }

  /// Get interpolated value for the first dimension of S(x). (Simplified interface for 1D)
  GPUd() DataT interpolate(DataT x1, DataT x2) const
  {
    return interpolateMath(mYdim, x1, x2);
  }

  /// _______________  IO   ________________________

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE)
  /// write a class object to the file
  int writeToFile(TFile& outf, const char* name);

  /// read a class object from the file
  static Spline2DBase* readFromFile(TFile& inpf, const char* name);
#endif

  /// _______________  Getters   ________________________

  /// Get number of Y dimensions
  GPUd() int getYdimensions() const { return mYdim; }

  /// Get minimal required alignment for the spline parameters
  GPUd() static constexpr size_t getParameterAlignmentBytes() { return 16; }

  /// Number of parameters
  GPUd() int getNumberOfParameters() const { return this->getNumberOfParametersMath(mYdim); }

  /// Size of the parameter array in bytes
  GPUd() size_t getSizeOfParameters() const { return sizeof(DataT) * this->getNumberOfParameters(); }

  /// Get a number of knots
  GPUd() int getNumberOfKnots() const { return mGridU1.getNumberOfKnots() * mGridU2.getNumberOfKnots(); }

  /// Get 1-D grid for U1 coordinate
  GPUd() const Spline1D<DataT>& getGridU1() const { return mGridU1; }

  /// Get 1-D grid for U2 coordinate
  GPUd() const Spline1D<DataT>& getGridU2() const { return mGridU2; }

  /// Get 1-D grid for U1 or U2 coordinate
  GPUd() const Spline1D<DataT>& getGrid(int iu) const { return (iu == 0) ? mGridU1 : mGridU2; }

  /// Get u of i-th knot
  GPUd() void getKnotU(int iKnot, int& u1, int& u2) const;

  /// Get index of a knot (iKnotU1,iKnotU2)
  GPUd() int getKnotIndex(int iKnotU1, int iKnotU2) const;

  /// Get spline parameters
  GPUd() DataT* getParameters() { return mParameters; }

  /// Get spline parameters const
  GPUd() const DataT* getParameters() const { return mParameters; }

  /// _______________  Technical stuff  ________________________

  /// Get offset of GridU flat data in the flat buffer
  GPUd() size_t getGridU1Offset() const { return mGridU1.getFlatBufferPtr() - mFlatBufferPtr; }

  /// Get offset of GridU2 flat data in the flat buffer
  GPUd() size_t getGridU2Offset() const { return mGridU2.getFlatBufferPtr() - mFlatBufferPtr; }

  /// Set X range
  GPUd() void setXrange(DataT x1Min, DataT x1Max, DataT x2Min, DataT x2Max);

  /// Print method
  void print() const;

  ///  _______________  Expert tools  _______________

  /// Get interpolated value for an mYdim-dimensional S(u) using spline parameters Parameters.
  GPUd() void interpolateU(GPUgeneric() const DataT Parameters[],
                           DataT u1, DataT u2, GPUgeneric() DataT S[/*mYdim*/]) const
  {
    interpolateUMath(mYdim, Parameters, u1, u2, S);
  }

  /// Same as interpolateU(..) but without boundary checks
  GPUd() void interpolateUnonSafe(GPUgeneric() const DataT Parameters[],
                                  DataT u1, DataT u2, GPUgeneric() DataT S[/*mYdim*/]) const
  {
    interpolateUnonSafeMath(mYdim, Parameters, u1, u2, S);
  }

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE) // code invisible on GPU and in the standalone compilation
  /// Test the class functionality
  static int test(const bool draw = 0, const bool drawDataPoints = 1);
#endif

  /// _____________  FlatObject functionality, see FlatObject class for description  ____________

  using FlatObject::getBufferAlignmentBytes;
  using FlatObject::getClassAlignmentBytes;

#if !defined(GPUCA_GPUCODE)
  void cloneFromObject(const Spline2DBase& obj, char* newFlatBufferPtr);
  void moveBufferTo(char* newBufferPtr);
#endif

  using FlatObject::releaseInternalBuffer;

  void destroy();
  void setActualBufferAddress(char* actualFlatBufferPtr);
  void setFutureBufferAddress(char* futureFlatBufferPtr);

  /// _____________  Interpolation math with external Parameters ____________

  /// Number of parameters
  GPUd() int getNumberOfParametersMath(int nYdim) const { return (4 * nYdim) * getNumberOfKnots(); }

  /// Get interpolated value for an nYdim-dimensional S(u) using spline parameters Parameters.
  GPUd() void interpolateUMath(int nYdim, GPUgeneric() const DataT Parameters[],
                               DataT u1, DataT u2, GPUgeneric() DataT S[/*nYdim*/]) const;

  /// Same as interpolateU(..) but with no boundary checks
  GPUd() void interpolateUnonSafeMath(int nYdim, GPUgeneric() const DataT Parameters[],
                                      DataT u1, DataT u2, GPUgeneric() DataT S[/*nYdim*/]) const;

 protected:
  /// _____________  Interpolation math  ____________

  /// The main mathematical utility.
  /// Get interpolated value {S(u): 2D -> mYdim} at (u1,u2)
  /// using the spline values at the segment [knotL1, knotL1+1][knotL2, knotL2+1]
  /// Only first nYdimCalc out of nYdim Y dimensions are calculated
  template <SafeFlag Safe, bool DoFirstDimOnly>
  GPUd() void interpolateUMath1(int nYdim, GPUgeneric() const DataT Parameters[],
                                DataT u1, DataT u2, GPUgeneric() DataT S[/*nYdim*/]) const;

  /// Get interpolated value S(x)
  GPUd() void interpolateMath(int nYdim, DataT x1, DataT x2, GPUgeneric() DataT S[/*nYdim*/]) const;

  /// Get interpolated value for the first dimension of S(x). (Simplified interface for 1D)
  GPUd() DataT interpolateMath(int nYdim, DataT x1, DataT x2) const;

  /// _____________  Data members  ____________

 private:
  int mYdim; ///< dimentionality of F

 protected:                /// _____________  Data members  ____________
  Spline1D<DataT> mGridU1; ///< grid for U axis
  Spline1D<DataT> mGridU2; ///< grid for V axis
  DataT* mParameters;      //! (transient!!) F-dependent parameters of the spline
#ifndef GPUCA_ALIROOT_LIB
  ClassDefNV(Spline2DBase, 1);
#endif
};

///
/// An additional class specification where the spline dimensionality is set at the compile time
/// via a template argument nYdimT
/// It inherits from the main specification above and only replaces some methods by their faster versions.
///
template <class DataT, int nYdimT>
class Spline2D : public Spline2DBase<DataT>
{
 public:
  typedef Spline2DBase<DataT> TBase;

  /// _____________  C++ constructors / destructors __________________________

#if !defined(GPUCA_GPUCODE)
  /// Default constructor
  Spline2D() : Spline2D(2, 2) {}

  /// Constructor for a regular spline && default constructor
  Spline2D(int numberOfKnotsU1, int numberOfKnotsU2) : TBase(nYdimT, numberOfKnotsU1, numberOfKnotsU2) {}

  /// Constructor for an irregular spline
  Spline2D(int numberOfKnotsU1, const int knotU1[], int numberOfKnotsU2, const int knotU2[])
    : TBase(nYdimT, numberOfKnotsU1, knotU1, numberOfKnotsU2, knotU2) {}

  /// Copy constructor
  Spline2D(const Spline2D& v) : TBase(v) {}

  /// Assignment operator
  Spline2D& operator=(const Spline2D& v) { return (Spline2D&) TBase::operator=(v); };
#else
  /// Disable constructors for the GPU implementation
  Spline2D() CON_DELETE;
  Spline2D(const Spline2D&) CON_DELETE;
  Spline2D& operator=(const Spline2D&) CON_DELETE;
#endif

  /// Destructor
  ~Spline2D() CON_DEFAULT;

  /// _____________  Construction interface __________________________

#if !defined(GPUCA_GPUCODE)
  /// Constructor for a regular spline.
  void recreate(int numberOfKnotsU1, int numberOfKnotsU2) { TBase::recreate(nYdimT, numberOfKnotsU1, numberOfKnotsU2); }

  /// Constructor for an irregular spline
  void recreate(int numberOfKnotsU1, const int knotU1[], int numberOfKnotsU2, const int knotU2[])
  {
    TBase::recreate(nYdimT, numberOfKnotsU1, knotU1, numberOfKnotsU2, knotU2);
  }
#endif

  /// _______________  Main functionality   ________________________

  /// Get interpolated value S(x)
  GPUd() void interpolate(DataT x1, DataT x2, GPUgeneric() DataT S[/*mYdim*/]) const
  {
    interpolateMath(nYdimT, x1, x2, S);
  }

  /// Same as interpolate(), but using vectorized calculations
  GPUd() void interpolateVec(DataT x1, DataT x2, GPUgeneric() DataT S[/*mYdim*/]) const
  {
    interpolateMath(nYdimT, x1, x2, S);
  }

  /// Get interpolated value for the first dimension of S(x). (Simplified interface for 1D)
  GPUd() DataT interpolate(DataT x1, DataT x2) const
  {
    return interpolateMath(nYdimT, x1, x2);
  }

  /// _______________  Getters   ________________________

  /// Get number of F dimensions as constexpr
  GPUd() static constexpr int getYdimensions() { return nYdimT; }

  /// Number of parameters
  GPUd() int getNumberOfParameters() const { return this->getNumberOfParametersMath(nYdimT); }

  /// Size of the parameter array in bytes
  GPUd() size_t getSizeOfParameters() const { return sizeof(DataT) * this->getNumberOfParameters(); }

  ///  _______________  Expert tools  _______________

  /// Get interpolated value for an mYdim-dimensional S(u) using spline parameters Parameters.
  GPUd() void interpolateU(GPUgeneric() const DataT Parameters[],
                           DataT u1, DataT u2, GPUgeneric() DataT S[/*mYdim*/]) const
  {
    interpolateUMath(nYdimT, Parameters, u1, u2, S);
  }

  /// Same as interpolateU(..) but without boundary checks
  GPUd() void interpolateUnonSafe(GPUgeneric() const DataT Parameters[],
                                  DataT u1, DataT u2, GPUgeneric() DataT S[/*mYdim*/]) const
  {
    interpolateUnonSafeMath(nYdimT, Parameters, u1, u2, S);
  }

  /// _______________  IO   ________________________

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE)
  /// read a class object from the file
  static Spline2D* readFromFile(TFile& inpf, const char* name)
  {
    return reinterpret_cast<Spline2D*>(TBase::readFromFile(inpf, name));
  }
#endif

  /// _______________  Suppress some base class methods   ________________________
 private:
#if !defined(GPUCA_GPUCODE)
  void recreate(int nYdim, int numberOfKnotsU1, int numberOfKnotsU2){};
  void recreate(int nYdim, int numberOfKnotsU1, const int knotU1[], int numberOfKnotsU2, const int knotU2[]){};
#endif
 public:
  using TBase::getNumberOfKnots;
  using TBase::interpolateMath;
  using TBase::interpolateUMath;
  using TBase::interpolateUnonSafeMath;
  using TBase::mGridU1;
  using TBase::mGridU2;
  using TBase::mParameters;
};

///
/// ========================================================================================================
///       Inline implementations of some methods
/// ========================================================================================================
///

#if !defined(GPUCA_GPUCODE)
template <typename DataT>
GPUdi() Spline2DBase<DataT>::Spline2DBase(int nYdim,
                                          int numberOfKnotsU1, int numberOfKnotsU2)
  : FlatObject(),
    mYdim(0),
    mGridU1(),
    mGridU2(),
    mParameters(nullptr)
{ /// Constructor for a regular spline
  recreate(nYdim, numberOfKnotsU1, numberOfKnotsU2);
}

template <typename DataT>
GPUdi() Spline2DBase<DataT>::Spline2DBase(int nYdim,
                                          int numberOfKnotsU1, const int knotU1[],
                                          int numberOfKnotsU2, const int knotU2[])
  : FlatObject(),
    mYdim(0),
    mGridU1(),
    mGridU2(),
    mParameters(nullptr)
{ /// Constructor for an irregular spline
  recreate(nYdim, numberOfKnotsU1, knotU1, numberOfKnotsU2, knotU2);
}

template <typename DataT>
GPUdi() Spline2DBase<DataT>::Spline2DBase(const Spline2DBase& spline)
  : FlatObject(),
    mYdim(0),
    mGridU1(),
    mGridU2(),
    mParameters(nullptr)
{ /// Copy constructor
  cloneFromObject(spline, nullptr);
}

template <typename DataT>
GPUdi() Spline2DBase<DataT>& Spline2DBase<DataT>::operator=(const Spline2DBase<DataT>& spline)
{ /// Assignment operator
  cloneFromObject(spline, nullptr);
  return *this;
}
#endif

template <typename DataT>
GPUdi() void Spline2DBase<DataT>::getKnotU(int iKnot, int& u1, int& u2) const
{ /// Get u1,u2 of i-th knot
  int nu1 = mGridU1.getNumberOfKnots();
  int iu2 = iKnot / nu1;
  int iu1 = iKnot % nu1;
  u1 = mGridU1.getKnot(iu1).getU();
  u2 = mGridU2.getKnot(iu2).getU();
}

template <typename DataT>
GPUdi() int Spline2DBase<DataT>::getKnotIndex(int iKnotU1, int iKnotU2) const
{ /// Get index of a knot (iKnotU1,iKnotU2)
  int nu1 = mGridU1.getNumberOfKnots();
  return nu1 * iKnotU2 + iKnotU1;
}

template <typename DataT>
GPUdi() void Spline2DBase<DataT>::setXrange(
  DataT x1Min, DataT x1Max, DataT x2Min, DataT x2Max)
{ /// Set X range
  mGridU1.setXrange(x1Min, x1Max);
  mGridU2.setXrange(x2Min, x2Max);
}

template <typename DataT>
GPUdi() void Spline2DBase<DataT>::interpolateUMath(int nYdim, GPUgeneric() const DataT Parameters[],
                                                   DataT u1, DataT u2, GPUgeneric() DataT S[/*nYdim*/]) const
{
  /// Get interpolated value for an nYdim-dimensional S(u) using spline parameters Parameters.
  interpolateUMath1<SafeFlag::kSafe, 0>(nYdim, Parameters, u1, u2, S);
}

template <typename DataT>
GPUdi() void Spline2DBase<DataT>::interpolateUnonSafeMath(int nYdim, GPUgeneric() const DataT Parameters[],
                                                          DataT u1, DataT u2, GPUgeneric() DataT S[/*nYdim*/]) const
{
  /// Same as interpolateU(..) but with no boundary checks
  interpolateUMath1<SafeFlag::kNotSafe, 0>(nYdim, Parameters, u1, u2, S);
}

/// The main mathematical utility.
/// Get interpolated value {S(u): 2D -> mYdim} at (u1,u2)
/// using the spline values at the segment [knotL1, knotL1+1][knotL2, knotL2+1]
/// Only first nYdimCalc out of nYdim Y dimensions are calculated
template <typename DataT>
template <typename Spline2DBase<DataT>::SafeFlag Safe, bool DoFirstDimOnly>
GPUdi() void Spline2DBase<DataT>::interpolateUMath1(int nYdim, GPUgeneric() const DataT Parameters[],
                                                    DataT u1, DataT u2, GPUgeneric() DataT S[/*nYdim*/]) const
{
  /*
  float u = u1;
  float v = u2;
  int nu = mGridU1.getNumberOfKnots();
  int iu = mGridU1.getLeftKnotIndexForU<Safe>(u);
  int iv = mGridU2.getLeftKnotIndexForU<Safe>(v);

  const typename Spline1D<DataT>::Knot& knotU = mGridU1.getKnot<SafeFlag::kNotSafe>(iu);
  const typename Spline1D<DataT>::Knot& knotV = mGridU2.getKnot<SafeFlag::kNotSafe>(iv);

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

template <typename DataT>
GPUdi() void Spline2DBase<DataT>::
  interpolateMath(int nYdim, DataT x1, DataT x2, GPUgeneric() DataT S[/*nYdim*/]) const
{ /// Get interpolated value S(x)
  interpolateUMath1<SafeFlag::kSafe, 0>(nYdim, mParameters, mGridU1.convXtoU(x1), mGridU2.convXtoU(x2), S);
}

template <typename DataT>
GPUdi() DataT Spline2DBase<DataT>::
  interpolateMath(int nYdim, DataT x1, DataT x2) const
{ /// Get interpolated value for the first dimension of S(x). (Simplified interface for 1D)
  DataT S;
  interpolateUMath1<SafeFlag::kSafe, 1>(nYdim, mParameters, mGridU1.convXtoU(x1), mGridU2.convXtoU(x2), &S);
  return S;
}

#ifdef XXXX

template <typename DataT, int nFdimT>
GPUhdi() void Spline2D<DataT, nFdimT>::interpolateVec(
  DataT x1, DataT x2, GPUgeneric() DataT S[]) const
{
  /// Same as interpolate(), but using vectorized calculation
  interpolateUvec(mFparameters, mGridU1.convXtoU(x1), mGridU2.convXtoU(x2), S);
}

template <typename DataT, int nFdimT>
GPUhdi() void Spline2D<DataT, nFdimT>::interpolateU(
  GPUgeneric() const DataT Fparameters[],
  DataT u, DataT v, GPUgeneric() DataT S[]) const
{
  /// Get interpolated value for an nFdim-dimensional F(u) using spline parameters Fparameters.
  /// Fparameters can be created via SplineHelper2D.

  int nu = mGridU1.getNumberOfKnots();
  int iu = mGridU1.getLeftKnotIndexForU<>(u);
  int iv = mGridU2.getLeftKnotIndexForU<>(v);

  const typename Spline1D<DataT>::Knot& knotU = mGridU1.getKnot(iu);
  const typename Spline1D<DataT>::Knot& knotV = mGridU2.getKnot(iv);

#if defined(GPUCA_GPUCODE) // constexpr array size for the GPU compiler
  const int nFdim = nFdimT;
#else
  const int nFdim = getFdimensions();
#endif

  const int nFdim2 = nFdimT * 2;
  const int nFdim4 = nFdimT * 4;

  // X:=Sx, Y:=Sy, Z:=Sz

  const DataT* par00 = Fparameters + (nu * iv + iu) * nFdim4; // values { {X,Y,Z}, {X,Y,Z}'v, {X,Y,Z}'u, {X,Y,Z}''vu } at {u0, v0}
  const DataT* par10 = par00 + nFdim4;                        // values { ... } at {u1, v0}
  const DataT* par01 = par00 + nFdim4 * nu;                   // values { ... } at {u0, v1}
  const DataT* par11 = par01 + nFdim4;                        // values { ... } at {u1, v1}

  DataT Su0[nFdim4]; // values { {X,Y,Z,X'v,Y'v,Z'v}(v0), {X,Y,Z,X'v,Y'v,Z'v}(v1) }, at u0
  DataT Du0[nFdim4]; // derivatives {}'_u  at u0
  DataT Su1[nFdim4]; // values { {X,Y,Z,X'v,Y'v,Z'v}(v0), {X,Y,Z,X'v,Y'v,Z'v}(v1) }, at u1
  DataT Du1[nFdim4]; // derivatives {}'_u  at u1

  for (int i = 0; i < nFdim2; i++) {
    Su0[i] = par00[i];
    Su0[nFdim2 + i] = par01[i];

    Du0[i] = par00[nFdim2 + i];
    Du0[nFdim2 + i] = par01[nFdim2 + i];

    Su1[i] = par10[i];
    Su1[nFdim2 + i] = par11[i];

    Du1[i] = par10[nFdim2 + i];
    Du1[nFdim2 + i] = par11[nFdim2 + i];
  }

  DataT parU[nFdim4]; // interpolated values { {X,Y,Z,X'v,Y'v,Z'v}(v0), {X,Y,Z,X'v,Y'v,Z'v}(v1) } at u
  mGridU1.interpolateUMath(nFdim4, knotU, Su0, Du0, Su1, Du1, u, parU);

  const DataT* Sv0 = parU + 0;
  const DataT* Dv0 = parU + nFdim;
  const DataT* Sv1 = parU + nFdim2;
  const DataT* Dv1 = parU + nFdim2 + nFdim;

  mGridU2.interpolateUMath(nFdim, knotV, Sv0, Dv0, Sv1, Dv1, v, S);
}

template <typename DataT, int nFdimT>
GPUhdi() void Spline2D<DataT, nFdimT>::interpolateUvec(
  GPUgeneric() const DataT Fparameters[],
  DataT u1, DataT u2, GPUgeneric() DataT S[]) const
{
  /// Same as interpolateU(), but using vectorized calculation
  interpolateU(mFparameters, u1, u2, S);
}

#endif

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
