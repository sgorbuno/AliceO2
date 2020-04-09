// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  Spline.h
/// \brief Definition of Spline class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#ifndef ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SPLINE_H
#define ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SPLINE_H

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
/// The Spline class performs a cubic spline interpolation on a multi-dimensional nonunifom grid.
/// The class is an extension of the Spline1D class.
/// See Spline1D.h for more details.
///
/// The spline S(x) approximates a function F(x):R^n->R^m,
/// where x is a vector which belongs to some [xMin,xMax] n-dimensional area. F value may be also multi-dimensional.
///
/// --- Example of creating a spline ---
///
///  constexpr int nXdim=2, nFdim=1;
///  auto F = [&](float[nXdim] x, float f[nFdim] ) {
///   f[0] = 1.f + x[0] + x[1]*x[1]; // F(x)
///  };
///  const int nKnots[nXdim]={2,3};
///  int knotsU[nKnotsU] = {0, 1};
///  int knotsV[nKnotsV] = {0, 2, 5};
///  Spline<float,2,1> spline(nKnotsU, knotsU, nKnotsV, knotsV ); // prepare memory for 1-dimensional F
///  spline.approximateFunction(0., 1., 0.,1., F); //initialize spline to approximate F on area [0., 1.]x[0., 1.]
///  float S = spline.interpolate(.1, .3 ); // interpolated value at (.1,.3)
///
///  --- See also Spline::test();
///

/// Base class to store data members and non-inline methods
template <typename Tfloat, bool TisConsistent>
class SplineBase : public FlatObject
{
 public:
  /// _____________  Version control __________________________

  /// Version control
  GPUhd() static constexpr int getVersion() { return 1; }

  /// _____________  Constructors / destructors __________________________

#if !defined(GPUCA_GPUCODE)
  /// constructor && default constructor for the ROOT streamer
  SplineBase(int nDim = 1);
#else
  /// Disable constructors
  SplineBase() CON_DELETE;
#endif
  /// Disable constructors
  SplineBase(const SplineBase&) CON_DELETE;
  SplineBase& operator=(const SplineBase&) CON_DELETE;

  /// Destructor
  ~SplineBase() CON_DEFAULT;

  /// _______________  Construction interface  ________________________

#if !defined(GPUCA_GPUCODE)
  /// Constructor for a regular spline.
  void recreate(int numberOfKnotsU1, int numberOfKnotsU2);

  /// Constructor for an irregular spline
  void recreate(int numberOfKnotsU1, const int knotsU1[], int numberOfKnotsU2, const int knotsU2[]);

  /// approximate a function F with this spline.
  void approximateFunction(Tfloat x1Min, Tfloat x1Max, Tfloat x2Min, Tfloat x2Max,
                           std::function<void(Tfloat x1, Tfloat x2, Tfloat f[])> F,
                           int nAxiliaryDataPointsU1 = 4, int nAxiliaryDataPointsU2 = 4);
#endif

  /// _______________  IO   ________________________

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE)
  /// write a class object to the file
  int writeToFile(TFile& outf, const char* name);

  /// read a class object from the file
  static SplineBase* readFromFile(TFile& inpf, const char* name);
#endif

  /// _______________  Getters   ________________________

  /// Get number of F dimensions
  GPUhd() int getFdimensions() const { return mFdim; }

  ///
  GPUhd() static constexpr bool isConsistent() { return TisConsistent; }

  /// Get minimal required alignment for the spline parameters
  GPUhd() static constexpr size_t getParameterAlignmentBytes() { return 16; }

  /// Number of parameters
  GPUhd() int getNumberOfParameters() const { return (4 * mFdim) * getNumberOfKnots(); }

  /// Size of the parameter array in bytes
  GPUhd() size_t getSizeOfParameters() const { return sizeof(Tfloat) * getNumberOfParameters(); }

  /// Get number total of knots: UxV
  GPUhd() int getNumberOfKnots() const { return mGridU1.getNumberOfKnots() * mGridU2.getNumberOfKnots(); }

  /// Get 1-D grid for U1 coordinate
  GPUhd() const Spline1D<Tfloat>& getGridU1() const { return mGridU1; }

  /// Get 1-D grid for U2 coordinate
  GPUhd() const Spline1D<Tfloat>& getGridU2() const { return mGridU2; }

  /// Get 1-D grid for U1 or U2 coordinate
  GPUhd() const Spline1D<Tfloat>& getGrid(int iu) const { return (iu == 0) ? mGridU1 : mGridU2; }

  /// Get u1,u2 of i-th knot
  GPUhd() void getKnotU(int iKnot, Tfloat& u1, Tfloat& u2) const;

  /// Get index of a knot (iKnotU1,iKnotU2)
  GPUhd() int getKnotIndex(int iKnotU1, int iKnotU2) const;

  /// Get number of F parameters
  GPUhd() Tfloat* getFparameters() { return mFparameters; }

  /// Get number of F parameters
  GPUhd() const Tfloat* getFparameters() const { return mFparameters; }

  /// _______________  Technical stuff  ________________________

  /// Get offset of GridU flat data in the flat buffer
  size_t getGridU1Offset() const { return mGridU1.getFlatBufferPtr() - mFlatBufferPtr; }

  /// Get offset of GridU2 flat data in the flat buffer
  size_t getGridU2Offset() const { return mGridU2.getFlatBufferPtr() - mFlatBufferPtr; }

  /// Set X range
  void setXrange(Tfloat x1Min, Tfloat x1Max, Tfloat x2Min, Tfloat x2Max);

  /// Print method
  void print() const;

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE) // code invisible on GPU and in the standalone compilation
  /// Test the class functionality
  static int test(const bool draw = 0, const bool drawDataPoints = 1);
#endif

  /// _____________  FlatObject functionality, see FlatObject class for description  ____________

  using FlatObject::getBufferAlignmentBytes;
  using FlatObject::getClassAlignmentBytes;

#if !defined(GPUCA_GPUCODE)
  void cloneFromObject(const SplineBase& obj, char* newFlatBufferPtr);
  void moveBufferTo(char* newBufferPtr);
#endif

  using FlatObject::releaseInternalBuffer;

  void destroy();
  void setActualBufferAddress(char* actualFlatBufferPtr);
  void setFutureBufferAddress(char* futureFlatBufferPtr);

 private:
  int mFdim; ///< dimentionality of F

 protected:                 /// _____________  Data members  ____________
  Spline1D<Tfloat> mGridU1; ///< grid for U axis
  Spline1D<Tfloat> mGridU2; ///< grid for V axis
  Tfloat* mFparameters;     //! (transient!!) F-dependent parameters of the spline

  ClassDefNV(SplineBase, 1);
};

///
/// The main Spline class. Contains constructors and interpolation.
///
/// F dimensions can be set as a template parameter (faster; the only option for the GPU)
/// or as a constructor parameter (slower; the only option for the ROOT interpretator).
/// In a compiled CPU code one can use both options.
///
template <typename Tfloat, int TnFdim = 0, bool TisConsistent = 1>
class Spline : public SplineBase<Tfloat, TisConsistent>
{
 public:
  typedef SplineBase<Tfloat, TisConsistent> TBase;

  /// _____________  Constructors / destructors __________________________

#if !defined(GPUCA_GPUCODE)

  /// Constructor for a regular spline && default constructor
  Spline(int numberOfKnotsU1 = 2, int numberOfKnotsU2 = 2);

  /// Constructor for an irregular spline
  Spline(int numberOfKnotsU1, const int knotsU1[], int numberOfKnotsU2, const int knotsU2[]);

  /// Copy constructor
  Spline(const Spline&);

  /// Assignment operator
  Spline& operator=(const Spline&);
#else
  /// Disable constructors for the GPU implementation
  Spline() CON_DELETE;
  Spline(const Spline&) CON_DELETE;
  Spline& operator=(const Spline&) CON_DELETE;
#endif

  /// Destructor
  ~Spline() CON_DEFAULT;

  /// _______________  Main functionality   ________________________

  /// Get interpolated value for F(x1,x2)
  GPUhd() void interpolate(Tfloat x1, Tfloat x2, GPUgeneric() Tfloat S[]) const;

  /// Get interpolated value for the first dimension of F(x1,x2). (Simplified interface for 1D)
  GPUhd() Tfloat interpolate(Tfloat x1, Tfloat x2) const;

  /// Same as interpolate(), but using vectorized calculation.
  GPUhd() void interpolateVec(Tfloat x1, Tfloat x2, GPUgeneric() Tfloat S[]) const;

  /// ================ Expert tools   ================================

  /// Get interpolated value for an nFdim-dimensional F(u) using spline parameters Fparameters.
  /// Fparameters can be created via SplineHelper.
  GPUhd() void interpolateU(GPUgeneric() const Tfloat Fparameters[],
                            Tfloat u1, Tfloat u2, GPUgeneric() Tfloat Su[]) const;

  /// Same as interpolateU(), but using vectorized calculation.
  GPUhd() void interpolateUvec(GPUgeneric() const Tfloat Fparameters[],
                               Tfloat u1, Tfloat u2, GPUgeneric() Tfloat Su[]) const;

  /// _______________  IO   ________________________

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE)
  /// write a class object to the file
  using TBase::writeToFile;

  /// read a class object from the file
  static Spline* readFromFile(TFile& inpf, const char* name)
  {
    return reinterpret_cast<Spline*>(TBase::readFromFile(inpf, name));
  }
#endif

  /// _______________  Getters   ________________________

  /// Get number of F dimensions.
  GPUhd() static constexpr int getFdimensions() { return TnFdim; }

  using TBase::mFparameters;
  using TBase::mGridU1;
  using TBase::mGridU2;
};

///
/// ========================================================================================================
///       Inline implementations of some methods
/// ========================================================================================================
///

template <typename Tfloat, bool TisConsistent>
GPUhdi() void SplineBase<Tfloat, TisConsistent>::getKnotU(int iKnot, Tfloat& u1, Tfloat& u2) const
{
  /// Get u1,u2 of i-th knot
  int nu1 = mGridU1.getNumberOfKnots();
  int iu2 = iKnot / nu1;
  int iu1 = iKnot % nu1;
  u1 = mGridU1.getKnot(iu1).u;
  u2 = mGridU2.getKnot(iu2).u;
}

template <typename Tfloat, bool TisConsistent>
GPUhdi() int SplineBase<Tfloat, TisConsistent>::getKnotIndex(int iKnotU1, int iKnotU2) const
{
  /// Get index of a knot (iKnotU1,iKnotU2)
  int nu1 = mGridU1.getNumberOfKnots();
  return nu1 * iKnotU2 + iKnotU1;
}

template <typename Tfloat, bool TisConsistent>
GPUhdi() void SplineBase<Tfloat, TisConsistent>::setXrange(
  Tfloat x1Min, Tfloat x1Max, Tfloat x2Min, Tfloat x2Max)
{
  mGridU1.setXrange(x1Min, x1Max);
  mGridU2.setXrange(x2Min, x2Max);
}

#if !defined(GPUCA_GPUCODE)

template <typename Tfloat, int TnFdim, bool TisConsistent>
GPUhdi() Spline<Tfloat, TnFdim, TisConsistent>::
  Spline(int numberOfKnotsU1, int numberOfKnotsU2)
  : SplineBase<Tfloat, TisConsistent>(TnFdim)
{
  this->recreate(numberOfKnotsU1, numberOfKnotsU2);
}

template <typename Tfloat, int TnFdim, bool TisConsistent>
GPUhdi() Spline<Tfloat, TnFdim, TisConsistent>::
  Spline(int numberOfKnotsU1, const int knotsU1[],
         int numberOfKnotsU2, const int knotsU2[])
  : SplineBase<Tfloat, TisConsistent>(TnFdim)
{
  this->recreate(numberOfKnotsU1, knotsU1, numberOfKnotsU2, knotsU2);
}

template <typename Tfloat, int TnFdim, bool TisConsistent>
GPUhdi() Spline<Tfloat, TnFdim, TisConsistent>::
  Spline(const Spline& spline)
  : SplineBase<Tfloat, TisConsistent>(TnFdim)
{
  this->cloneFromObject(spline, nullptr);
}

template <typename Tfloat, int TnFdim, bool TisConsistent>
GPUhdi() Spline<Tfloat, TnFdim, TisConsistent>& Spline<Tfloat, TnFdim, TisConsistent>::
  operator=(const Spline& spline)
{
  this->cloneFromObject(spline, nullptr);
  return *this;
}
#endif

template <typename Tfloat, int TnFdim, bool TisConsistent>
GPUhdi() void Spline<Tfloat, TnFdim, TisConsistent>::
  interpolate(Tfloat x1, Tfloat x2, GPUgeneric() Tfloat S[]) const
{
  /// Get interpolated value for F(x1,x2)
  assert(TisConsistent);
  if (TisConsistent) {
    interpolateU(mFparameters, mGridU1.convXtoU(x1), mGridU2.convXtoU(x2), S);
  } else {
    for (int i = 0; i < TnFdim; i++) {
      S[i] = -1.;
    }
  }
}

template <typename Tfloat, int TnFdim, bool TisConsistent>
GPUhdi() Tfloat Spline<Tfloat, TnFdim, TisConsistent>::
  interpolate(Tfloat x1, Tfloat x2) const
{
  /// Simplified interface for 1D: get interpolated value for the first dimension of F(x)

#if defined(GPUCA_GPUCODE)
  Tfloat S[TnFdim]; // constexpr array size for the GPU compiler
#else
  Tfloat S[getFdimensions()];
#endif
  interpolate(x1, x2, S);
  return S[0];
}

template <typename Tfloat, int TnFdim, bool TisConsistent>
GPUhdi() void Spline<Tfloat, TnFdim, TisConsistent>::interpolateVec(
  Tfloat x1, Tfloat x2, GPUgeneric() Tfloat S[]) const
{
  /// Same as interpolate(), but using vectorized calculation
  assert(TisConsistent);
  if (TisConsistent) {
    interpolateUvec(mFparameters, mGridU1.convXtoU(x1), mGridU2.convXtoU(x2), S);
  } else {
    for (int i = 0; i < getFdimensions(); i++) {
      S[i] = 0.;
    }
  }
}

template <typename Tfloat, int TnFdim, bool TisConsistent>
GPUhdi() void Spline<Tfloat, TnFdim, TisConsistent>::interpolateU(
  GPUgeneric() const Tfloat Fparameters[],
  Tfloat u, Tfloat v, GPUgeneric() Tfloat S[]) const
{
  /// Get interpolated value for an nFdim-dimensional F(u) using spline parameters Fparameters.
  /// Fparameters can be created via SplineHelper.

  int nu = mGridU1.getNumberOfKnots();
  int iu = mGridU1.getKnotIndexU(u);
  int iv = mGridU2.getKnotIndexU(v);

  const typename Spline1D<Tfloat>::Knot& knotU = mGridU1.getKnot(iu);
  const typename Spline1D<Tfloat>::Knot& knotV = mGridU2.getKnot(iv);

#if defined(GPUCA_GPUCODE) // constexpr array size for the GPU compiler
  const int nFdim = TnFdim;
#else
  const int nFdim = getFdimensions();
#endif

  const int nFdim2 = TnFdim * 2;
  const int nFdim4 = TnFdim * 4;

  // X:=Sx, Y:=Sy, Z:=Sz

  const Tfloat* par00 = Fparameters + (nu * iv + iu) * nFdim4; // values { {X,Y,Z}, {X,Y,Z}'v, {X,Y,Z}'u, {X,Y,Z}''vu } at {u0, v0}
  const Tfloat* par10 = par00 + nFdim4;                        // values { ... } at {u1, v0}
  const Tfloat* par01 = par00 + nFdim4 * nu;                   // values { ... } at {u0, v1}
  const Tfloat* par11 = par01 + nFdim4;                        // values { ... } at {u1, v1}

  Tfloat Su0[nFdim4]; // values { {X,Y,Z,X'v,Y'v,Z'v}(v0), {X,Y,Z,X'v,Y'v,Z'v}(v1) }, at u0
  Tfloat Du0[nFdim4]; // derivatives {}'_u  at u0
  Tfloat Su1[nFdim4]; // values { {X,Y,Z,X'v,Y'v,Z'v}(v0), {X,Y,Z,X'v,Y'v,Z'v}(v1) }, at u1
  Tfloat Du1[nFdim4]; // derivatives {}'_u  at u1

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

  Tfloat parU[nFdim4]; // interpolated values { {X,Y,Z,X'v,Y'v,Z'v}(v0), {X,Y,Z,X'v,Y'v,Z'v}(v1) } at u
  mGridU1.interpolateU(nFdim4, knotU, Su0, Du0, Su1, Du1, u, parU);

  const Tfloat* Sv0 = parU + 0;
  const Tfloat* Dv0 = parU + nFdim;
  const Tfloat* Sv1 = parU + nFdim2;
  const Tfloat* Dv1 = parU + nFdim2 + nFdim;

  mGridU2.interpolateU(nFdim, knotV, Sv0, Dv0, Sv1, Dv1, v, S);
}

template <typename Tfloat, int TnFdim, bool TisConsistent>
GPUhdi() void Spline<Tfloat, TnFdim, TisConsistent>::interpolateUvec(
  GPUgeneric() const Tfloat Fparameters[],
  Tfloat u1, Tfloat u2, GPUgeneric() Tfloat S[]) const
{
  /// Same as interpolateU(), but using vectorized calculation
  interpolateU(mFparameters, u1, u2, S);
}

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
