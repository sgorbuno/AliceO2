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

/*Abdel*/ #include <fstream>
#include "Spline1D.h"
#include "FlatObject.h"
#include "GPUCommonDef.h"
#if !defined(__CINT__) && !defined(__ROOTCINT__) && !defined(GPUCA_GPUCODE) && !defined(GPUCA_NO_VC) && defined(__cplusplus) && __cplusplus >= 201703L
#include <Vc/Vc>
#include <Vc/SimdArray>
#endif

namespace GPUCA_NAMESPACE
{
namespace gpu
{
///
/// The Spline2D class represents spline interpolation on a two-dimensional nonunifom (irregular) grid.
/// The class is an extension of the Spline1D class, see Spline1D.h for more details.
///
/// The class is a flat C structure. No virtual methods, no ROOT types are used.
/// It is designed for the parameterisation of TPC transformation.
///
/// --- Interpolation ---
///
/// The spline interpolates a function F:[u,v]->R^m,
/// where u,v belong to [0,Umax]x[0,Vmax].
///
/// The interpolation is performed on knots {U0==0., U1, .., Un-1==Umax} x {V0==0., V1, .., Vn-1==Vmax}
/// using the function values Fi and the derivatives Fi'_u, Fi'v, Fi''vu at the knots.
///
/// The interpolation of F values between the knots is performed by the 3-th degree polynoms.
/// The polynoms and they 1-st derivatives are continuous at the knots.
/// Depending on the initialization of the derivatives, the second derivative may or may not be continious.
///
/// --- Knots ---
///
/// Umax, Vmax are integer numbers.
/// The knots belong to [0, Umax][0,Vmax] and have integer coordinates.
/// It is implemented this way for fast matching of UV coordinates to the corresponding U&V segments between the knots.
///
/// To interpolate a function on an interval other than [0,Umax]x[0,Vmax], one should scale the U/V coordinates and the F derivatives in the parameter array.
///
/// The minimal number of knots is 2, the minimal Umax & Vmax is 1
///
/// The knot U0=0 must always be present. It will be added automatically when it is not set by the user.
/// Duplicated knots are not allowed and will be removed.
/// Therefore: Attention! The number of knots may change during the initialization.
///
/// ---- The parameter array ---
///
/// Nothing which depends on F is stored in the class,
/// therefore one can use the same class for interpolation of different input functions on the same knots.
///
/// The user should create an array of the function values Fi and all the derivatives at all (initialised!) knots.
/// This parameter array can be created using utilities from the SplineHelper2D class.
///
///  The parameter array is an array with info about the interpolated function F:[u,v]->[Fx,Fy,Fz] at knots:
///   {
///     { (Fx,Fy,Fz), (Fx,Fy,Fz)'_v, (Fx,Fy,Fz)'_u, (Fx,Fy,Fz)''_vu } at knot 0,
///     {                      ...                      } at knot 1,
///                            ...
///   }
///   The parameter array has to be provided by the user for each call of the interpolation.
///   It can be created for a given input function using SplineHelper2D class.
///
/// ---- Flat Object implementation ----
///
/// The class inherits from the FlatObjects. Copying can only be done using the FlatObject interface.
///
/// --- Example of creating a spline ---
///
///  auto F = [&](float u, float v ) -> float {
///   return ...; // F(u,v)
///  };
///  const int nKnotsU=2;
///  const int nKnotsV=3;
///  int knotsU[nKnotsU] = {0, 2};
///  int knotsV[nKnotsV] = {0, 3, 6};
///  Spline2D spline;
///  spline.construct(nKnotsU, knotsU, nKnotsV, knotsV );
///  SplineHelper2D helper;
///  std::unique_ptr<float[]> parameters = helper.constructParameters(1,spline, F, 0.f, 1.f, 0.f, 1.f, 2);
///  spline.interpolate(1, parameters.get(), 0.0, 0.0 ); // == F(0.,0.)
///  spline.interpolate(1, parameters.get(), 1.0, 1.1); // some interpolated value
///  spline.interpolate(1, parameters.get(), 2.0, 3.0 ); // == F(1., 0.5 )
///  spline.interpolate(1, parameters.get(), 2.0, 6.0 ); // == F(1., 1.)
///
///  --- See also Spline2D::test();
///
class Spline2D : public FlatObject
{
 public:
  /// _____________  Version control __________________________

  /// Version number
  GPUd() static constexpr int getVersion() { return 1; }

  /// _____________  Constructors / destructors __________________________

  /// Default constructor. Creates an empty uninitialised object
  Spline2D();

  /// Copy constructor: disabled to avoid ambiguity. Use cloneFromObject() instead
  Spline2D(const Spline2D&) CON_DELETE;

  /// Assignment operator: disabled to avoid ambiguity. Use cloneFromObject() instead
  Spline2D& operator=(const Spline2D&) CON_DELETE;

  /// Destructor
  ~Spline2D() CON_DEFAULT;

  /// _____________  FlatObject functionality, see FlatObject class for description  ____________

  /// Memory alignment

  using FlatObject::getBufferAlignmentBytes;
  using FlatObject::getClassAlignmentBytes;

  /// Construction interface

  void cloneFromObject(const Spline2D& obj, char* newFlatBufferPtr);
  void destroy();

  /// Making the data buffer external

  using FlatObject::releaseInternalBuffer;
  void moveBufferTo(char* newBufferPtr);

  /// Moving the class with its external buffer to another location

  void setActualBufferAddress(char* actualFlatBufferPtr);
  void setFutureBufferAddress(char* futureFlatBufferPtr);

  /// _______________  Construction interface  ________________________

#if !defined(GPUCA_GPUCODE)
  /// Constructor
  ///
  /// Number of created knots may differ from the input values:
  /// - Edge knots {0} and {Umax/Vmax} will be added if they are not present.
  /// - Duplicated knots, knots with a negative coordinate will be deleted
  /// - At least 2 knots for each axis will be created
  ///
  /// \param numberOfKnotsU     Number of knots in knotsU[] array
  /// \param knotsU             Array of knot positions (integer values)
  ///
  /// \param numberOfKnotsV     Number of knots in knotsV[] array
  /// \param knotsV             Array of knot positions (integer values)
  ///
  void construct(int numberOfKnotsU, const int knotsU[], int numberOfKnotsV, const int knotsV[]);

  /// Constructor for a regular spline. Knots will be placed at the positions i/(numberOfKnots-1)
  void constructRegular(int numberOfKnotsU, int numberOfKnotsV);
#endif

  /// _______________  Main functionality   ________________________

  /// Get interpolated value for f(u,v)
  template <typename T>
  GPUd() void interpolate(int Ndim, GPUgeneric() const T* parameters, float u, float v, GPUgeneric() T Suv[/*Ndim*/]) const;
  template <typename T>
  GPUd() void interpolateOptimizedScalar(int Ndim, GPUgeneric() const T* parameters, float u, float v, GPUgeneric() T Suv[/*Ndim*/]) const;
  /// Same as interpolate, but using vectorized calculation.
  /// \param parameters should be at least 128-bit aligned
  /*Abdel*/ GPUdi() void interpolateVecVerticalAoV(int Ndim, GPUgeneric() const float* parameters, float u, float v, GPUgeneric() float Suv[/*Ndim*/]) const;
  /*Abdel*/ GPUdi() void interpolateVecVerticalSingleVector(int Ndim, GPUgeneric() const float* parameters, float u, float v, GPUgeneric() float Suv[/*Ndim*/]) const;
  /*Abdel*/ GPUd() void interpolateVecVertical(int Ndim, GPUgeneric() const float* parameters, float u, float v, GPUgeneric() float Suv[/*Ndim*/]) const;
  /*Abdel*/ GPUd() void interpolateVecVerticalTest(int Ndim, GPUgeneric() const float* parameters, float u, float v, GPUgeneric() float Suv[/*Ndim*/]) const;
  
  /*Abdel*/ GPUdi() void interpolateVecHorizontal(int Ndim, GPUgeneric() const float* parameters, Vc::float_v u, Vc::float_v v, GPUgeneric() Vc::float_v Suv[/*Vc::float_v::size()*/]) const;
  /*Abdel*/ GPUdi() void interpolateVecHorizontalV2(int Ndim, GPUgeneric() const float* parameters, Vc::float_v u, Vc::float_v v, GPUgeneric() Vc::float_v Suv[/*Ndim*/]) const;
  /*Abdel*/ GPUdi() void interpolateVecHorizontalV3(int Ndim, GPUgeneric() const float* parameters, Vc::float_v u, Vc::float_v v, GPUgeneric() Vc::float_v Suv[/*Ndim*/]) const;

  /// _______________  Getters   ________________________

  /// Get minimal required alignment for the spline parameters

  template <typename T>
  static constexpr size_t getParameterAlignmentBytes(int Ndim)
  {
    return std::min<4 * sizeof(T) * Ndim, 16>;
  }

  /// Size of the parameter array in bytes
  template <typename T>
  GPUd() size_t getSizeOfParameters(int Ndim) const
  {
    return sizeof(T) * (size_t)getNumberOfParameters(Ndim);
  }

  /// Number of parameters
  GPUd() int getNumberOfParameters(int Ndim) const
  {
    return (4 * Ndim) * getNumberOfKnots();
  }

  /// Get number total of knots: UxV
  GPUd() int getNumberOfKnots() const { return mGridU.getNumberOfKnots() * mGridV.getNumberOfKnots(); }

  /// Get 1-D grid for U coordinate
  GPUd() const Spline1D& getGridU() const { return mGridU; }

  /// Get 1-D grid for V coordinate
  GPUd() const Spline1D& getGridV() const { return mGridV; }

  /// Get 1-D grid for U or V coordinate
  GPUd() const Spline1D& getGrid(int uv) const { return (uv == 0) ? mGridU : mGridV; }

  /// Get u,v of i-th knot
  GPUd() void getKnotUV(int iKnot, float& u, float& v) const;

  /// _______________  Technical stuff  ________________________

  /// Get offset of GridU flat data in the flat buffer
  size_t getGridUOffset() const { return mGridU.getFlatBufferPtr() - mFlatBufferPtr; }

  /// Get offset of GridV flat data in the flat buffer
  size_t getGridVOffset() const { return mGridV.getFlatBufferPtr() - mFlatBufferPtr; }

  /// Print method
  void print() const;

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE) // code invisible on GPU and in the standalone compilation
  /// Test the class functionality
  static int test(bool draw = 0);
#endif

 private:
  ///
  /// ====  Data members   ====
  ///

  Spline1D mGridU; ///< grid for U axis
  Spline1D mGridV; ///< grid for V axis
};

/// ====================================================
///       Inline implementations of some methods
/// ====================================================

GPUdi() void Spline2D::getKnotUV(int iKnot, float& u, float& v) const
{
  /// Get u,v of i-th knot
  const Spline1D& gridU = getGridU();
  const Spline1D& gridV = getGridV();
  int nu = gridU.getNumberOfKnots();
  int iv = iKnot / nu;
  int iu = iKnot % nu;
  u = gridU.getKnot(iu).u;
  v = gridV.getKnot(iv).u;
}

template <typename T>                 
GPUdi() void Spline2D::interpolate(int Ndim, GPUgeneric() const T* parameters, float u, float v, GPUgeneric() T Suv[/*Ndim*/]) const
{
  // Get interpolated value for f(u,v) using parameters[getNumberOfParameters()]
  const Spline1D& gridU = getGridU();
  const Spline1D& gridV = getGridV();
  int nu = gridU.getNumberOfKnots();
  int iu = gridU.getKnotIndex(u);
  int iv = gridV.getKnotIndex(v);

  const Spline1D::Knot& knotU = gridU.getKnot(iu);
  const Spline1D::Knot& knotV = gridV.getKnot(iv);

  const int Ndim2 = Ndim * 2; 
  const int Ndim4 = Ndim * 4; 

  // X:=Sx, Y:=Sy, Z:=Sz
  const T* par00 = parameters + (nu * iv + iu) * Ndim4; // values { {X,Y,Z}, {X,Y,Z}'v, {X,Y,Z}'u, {X,Y,Z}''vu } at {u0, v0}
  const T* par10 = par00 + Ndim4;                       // values { ... } at {u1, v0}
  const T* par01 = par00 + Ndim4 * nu;                  // values { ... } at {u0, v1}
  const T* par11 = par01 + Ndim4;                       // values { ... } at {u1, v1}

  T Su0[Ndim4]; // values { {X,Y,Z,X'v,Y'v,Z'v}(v0), {X,Y,Z,X'v,Y'v,Z'v}(v1) }, at u0
  T Du0[Ndim4]; // derivatives {}'_u  at u0
  T Su1[Ndim4]; // values { {X,Y,Z,X'v,Y'v,Z'v}(v0), {X,Y,Z,X'v,Y'v,Z'v}(v1) }, at u1
  T Du1[Ndim4]; // derivatives {}'_u  at u1

  for (int i = 0; i < Ndim2; i++) {
    Su0[i] = par00[i];
    Su0[Ndim2 + i] = par01[i];

    Du0[i] = par00[Ndim2 + i];
    Du0[Ndim2 + i] = par01[Ndim2 + i];

    Su1[i] = par10[i];
    Su1[Ndim2 + i] = par11[i];

    Du1[i] = par10[Ndim2 + i];
    Du1[Ndim2 + i] = par11[Ndim2 + i];
  }

  T parU[Ndim4]; // interpolated values { {X,Y,Z,X'v,Y'v,Z'v}(v0), {X,Y,Z,X'v,Y'v,Z'v}(v1) } at u

  gridU.interpolate<T>(Ndim4, knotU, Su0, Du0, Su1, Du1, u, parU);

  const T* Sv0 = parU + 0;
  const T* Dv0 = parU + Ndim;
  const T* Sv1 = parU + Ndim2;
  const T* Dv1 = parU + Ndim2 + Ndim;

  gridV.interpolate<T>(Ndim, knotV, Sv0, Dv0, Sv1, Dv1, v, Suv);
}

template <typename T>                 
GPUdi() void Spline2D::interpolateOptimizedScalar(int Ndim, GPUgeneric() const T* parameters, float u, float v, GPUgeneric() T Suv[/*Ndim*/]) const
{
  // Get interpolated value for f(u,v) using parameters[getNumberOfParameters()]
  const Spline1D& gridU = getGridU();
  const Spline1D& gridV = getGridV();
  int nu = gridU.getNumberOfKnots();
  int iu = gridU.getKnotIndex(u);
  int iv = gridV.getKnotIndex(v);

  const Spline1D::Knot& knotU = gridU.getKnot(iu);
  const Spline1D::Knot& knotV = gridV.getKnot(iv);

  const int Ndim2 = Ndim * 2; 
  const int Ndim4 = Ndim * 4; 

  // X:=Sx, Y:=Sy, Z:=Sz
  const T* par00 = parameters + (nu * iv + iu) * Ndim4; // values { {X,Y,Z}, {X,Y,Z}'v, {X,Y,Z}'u, {X,Y,Z}''vu } at {u0, v0}
  const T* par10 = par00 + Ndim4;                       // values { ... } at {u1, v0}
  const T* par01 = par00 + Ndim4 * nu;                  // values { ... } at {u0, v1}
  const T* par11 = par01 + Ndim4;                       // values { ... } at {u1, v1}

  T Su0[Ndim4]; // values { {X,Y,Z,X'v,Y'v,Z'v}(v0), {X,Y,Z,X'v,Y'v,Z'v}(v1) }, at u0
  T Du0[Ndim4]; // derivatives {}'_u  at u0
  T Su1[Ndim4]; // values { {X,Y,Z,X'v,Y'v,Z'v}(v0), {X,Y,Z,X'v,Y'v,Z'v}(v1) }, at u1
  T Du1[Ndim4]; // derivatives {}'_u  at u1

  memcpy(Su0, par00, Ndim2 * sizeof(T));
  memcpy(Su0+Ndim2, par01, Ndim2 * sizeof(T));

  memcpy(Du0, par00+Ndim2, Ndim2 * sizeof(T));
  memcpy(Du0+Ndim2, par01+Ndim2, Ndim2 * sizeof(T));

  memcpy(Su1, par10, Ndim2 * sizeof(T));
  memcpy(Su1+Ndim2, par11, Ndim2 * sizeof(T));

  memcpy(Du1, par10+Ndim2, Ndim2 * sizeof(T));
  memcpy(Du1+Ndim2, par11+Ndim2, Ndim2 * sizeof(T));

  T parU[Ndim4]; // interpolated values { {X,Y,Z,X'v,Y'v,Z'v}(v0), {X,Y,Z,X'v,Y'v,Z'v}(v1) } at u

  gridU.interpolate<T>(Ndim4, knotU, Su0, Du0, Su1, Du1, u, parU);
  
  const T* Sv0 = parU + 0;
  const T* Dv0 = parU + Ndim;
  const T* Sv1 = parU + Ndim2;
  const T* Dv1 = parU + Ndim2 + Ndim;

  gridV.interpolate<T>(Ndim, knotV, Sv0, Dv0, Sv1, Dv1, v, Suv);
}

GPUdi() void Spline2D::interpolateVecVerticalAoV(int Ndim, GPUgeneric() const float* parameters, float u, float v, GPUgeneric() float Suv[/*Ndim*/]) const
{
  const Spline1D& gridU = getGridU();
  const Spline1D& gridV = getGridV();
  int nu = gridU.getNumberOfKnots();
  int iu = gridU.getKnotIndex(u);
  int iv = gridV.getKnotIndex(v);

  const Spline1D::Knot& knotU = gridU.getKnot(iu);
  const Spline1D::Knot& knotV = gridV.getKnot(iv);

  const int Ndim2 = Ndim * 2; 
  const int Ndim4 = Ndim * 4; 

  // X:=Sx, Y:=Sy, Z:=Sz
  const float* par00 = parameters + (nu * iv + iu) * Ndim4; // values { {X,Y,Z}, {X,Y,Z}'v, {X,Y,Z}'u, {X,Y,Z}''vu } at {u0, v0}
  const float* par10 = par00 + Ndim4;                       // values { ... } at {u1, v0}
  const float* par01 = par00 + Ndim4 * nu;                  // values { ... } at {u0, v1}
  const float* par11 = par01 + Ndim4;                       // values { ... } at {u1, v1}

  Vc::float_v Su0[3], Du0[3], Su1[3], Du1[3], parU[3];  // Array Of Vector

  memcpy(Su0, par00, Ndim2 * sizeof(float));
  memcpy(((float*)(Su0)) + Ndim2, par01, Ndim2 * sizeof(float));  

  memcpy(Du0, par00 + Ndim2, Ndim2 * sizeof(float));
  memcpy(((float*)(Du0)) + Ndim2, par01 + Ndim2, Ndim2 * sizeof(float));

  memcpy(Su1, par10, Ndim2 * sizeof(float));
  memcpy(((float*)(Su1)) + Ndim2, par11, Ndim2 * sizeof(float));

  memcpy(Du1, par10 + Ndim2, Ndim2 * sizeof(float));
  memcpy(((float*)(Du1)) + Ndim2, par11 + Ndim2, Ndim2 * sizeof(float));

  gridU.interpolate(knotU, Su0[0], Du0[0], Su1[0], Du1[0], u, parU[0]);
  gridU.interpolate(knotU, Su0[1], Du0[1], Su1[1], Du1[1], u, parU[1]);
  gridU.interpolate(knotU, Su0[2], Du0[2], Su1[2], Du1[2], u, parU[2]);

  Vc::float_v Sv0(((float*)(parU)), Vc::Aligned), Dv0(((float*)(parU)) + 3), Sv1(((float*)(parU)) + 6), Dv1(((float*)(parU)) + 9), SuVec;
  gridV.interpolate(knotV, Sv0, Dv0, Sv1, Dv1, v, SuVec);
  Suv[0] = SuVec[0];
  Suv[1] = SuVec[1];
  Suv[2] = SuVec[2];
}

GPUdi() void Spline2D::interpolateVecVerticalSingleVector(int Ndim, GPUgeneric() const float* parameters, float u, float v, GPUgeneric() float Suv[/*Ndim*/]) const
{
  // Get interpolated value for f(u,v) using parameters[getNumberOfParameters()]
  const Spline1D& gridU = getGridU();
  const Spline1D& gridV = getGridV();
  int nu = gridU.getNumberOfKnots();
  int iu = gridU.getKnotIndex(u);
  int iv = gridV.getKnotIndex(v);

  const Spline1D::Knot& knotU = gridU.getKnot(iu);
  const Spline1D::Knot& knotV = gridV.getKnot(iv);

  const int Ndim2 = Ndim * 2; 
  const int Ndim4 = Ndim * 4; 

  // X:=Sx, Y:=Sy, Z:=Sz
  const float* par00 = parameters + (nu * iv + iu) * Ndim4; // values { {X,Y,Z}, {X,Y,Z}'v, {X,Y,Z}'u, {X,Y,Z}''vu } at {u0, v0}
  const float* par10 = par00 + Ndim4;                       // values { ... } at {u1, v0}
  const float* par01 = par00 + Ndim4 * nu;                  // values { ... } at {u0, v1}
  const float* par11 = par01 + Ndim4;                       // values { ... } at {u1, v1}

  float Su0[Ndim4]; // values { {X,Y,Z,X'v,Y'v,Z'v}(v0), {X,Y,Z,X'v,Y'v,Z'v}(v1) }, at u0
  float Du0[Ndim4]; // derivatives {}'_u  at u0
  float Su1[Ndim4]; // values { {X,Y,Z,X'v,Y'v,Z'v}(v0), {X,Y,Z,X'v,Y'v,Z'v}(v1) }, at u1
  float Du1[Ndim4]; // derivatives {}'_u  at u1

  memcpy(Su0, par00, Ndim2 * sizeof(float));
  memcpy(Su0+Ndim2, par01, Ndim2 * sizeof(float));

  memcpy(Du0, par00+Ndim2, Ndim2 * sizeof(float));
  memcpy(Du0+Ndim2, par01+Ndim2, Ndim2 * sizeof(float));

  memcpy(Su1, par10, Ndim2 * sizeof(float));
  memcpy(Su1+Ndim2, par11, Ndim2 * sizeof(float));

  memcpy(Du1, par10+Ndim2, Ndim2 * sizeof(float));
  memcpy(Du1+Ndim2, par11+Ndim2, Ndim2 * sizeof(float));

  Vc::float_v Su0Vec(Su0, Vc::Aligned), Du0Vec(Du0, Vc::Aligned), Su1Vec(Su1, Vc::Aligned), Du1Vec(Du1, Vc::Aligned), parUVec;
  float parU[12];

  gridU.interpolate(knotU, Su0Vec, Du0Vec, Su1Vec, Du1Vec, u, parUVec);
  parUVec.store(parU, Vc::Aligned);
  Su0Vec = Vc::float_v(Su0+4, Vc::Aligned); Du0Vec = Vc::float_v(Du0+4, Vc::Aligned); Su1Vec = Vc::float_v(Su1+4, Vc::Aligned); Du1Vec = Vc::float_v(Du1+4, Vc::Aligned);
  gridU.interpolate(knotU, Su0Vec, Du0Vec, Su1Vec, Du1Vec, u, parUVec);
  parUVec.store(parU+4, Vc::Aligned);
  Su0Vec = Vc::float_v(Su0+8, Vc::Aligned); Du0Vec = Vc::float_v(Du0+8, Vc::Aligned); Su1Vec = Vc::float_v(Su1+8, Vc::Aligned); Du1Vec = Vc::float_v(Du1+8, Vc::Aligned);
  gridU.interpolate(knotU, Su0Vec, Du0Vec, Su1Vec, Du1Vec, u, parUVec);
  parUVec.store(parU+8, Vc::Aligned);

  Vc::float_v Sv0Vec(parU, Vc::Aligned), Dv0Vec(parU + 3), Sv1Vec(parU + 6), Dv1Vec(parU + 9), SuVec;
  gridV.interpolate(knotV, Sv0Vec, Dv0Vec, Sv1Vec, Dv1Vec, v, SuVec);
  Suv[0] = SuVec[0];
  Suv[1] = SuVec[1];
  Suv[2] = SuVec[2];
}

GPUdi() void Spline2D::interpolateVecVerticalTest(int Ndim, GPUgeneric() const float* parameters, float u, float v, GPUgeneric() float Suv[/*Ndim*/]) const
{
  // Get interpolated value for f(u,v) using parameters[getNumberOfParameters()]
  const Spline1D& gridU = getGridU();
  const Spline1D& gridV = getGridV();
  int nu = gridU.getNumberOfKnots();
  int iu = gridU.getKnotIndex(u);
  int iv = gridV.getKnotIndex(v);

  const Spline1D::Knot& knotU = gridU.getKnot(iu);
  const Spline1D::Knot& knotV = gridV.getKnot(iv);

  const int Ndim2 = Ndim * 2; 
  const int Ndim4 = Ndim * 4; 

  // X:=Sx, Y:=Sy, Z:=Sz
  const float* par00 = parameters + (nu * iv + iu) * Ndim4; // values { {X,Y,Z}, {X,Y,Z}'v, {X,Y,Z}'u, {X,Y,Z}''vu } at {u0, v0}
  const float* par10 = par00 + Ndim4;                       // values { ... } at {u1, v0}
  const float* par01 = par00 + Ndim4 * nu;                  // values { ... } at {u0, v1}
  const float* par11 = par01 + Ndim4;                       // values { ... } at {u1, v1}

  float Su0[Ndim4]; // values { {X,Y,Z,X'v,Y'v,Z'v}(v0), {X,Y,Z,X'v,Y'v,Z'v}(v1) }, at u0
  float Du0[Ndim4]; // derivatives {}'_u  at u0
  float Su1[Ndim4]; // values { {X,Y,Z,X'v,Y'v,Z'v}(v0), {X,Y,Z,X'v,Y'v,Z'v}(v1) }, at u1
  float Du1[Ndim4]; // derivatives {}'_u  at u1

  memcpy(Su0, par00, Ndim2 * sizeof(float));
  memcpy(Su0+Ndim2, par01, Ndim2 * sizeof(float));

  memcpy(Du0, par00+Ndim2, Ndim2 * sizeof(float));
  memcpy(Du0+Ndim2, par01+Ndim2, Ndim2 * sizeof(float));

  memcpy(Su1, par10, Ndim2 * sizeof(float));
  memcpy(Su1+Ndim2, par11, Ndim2 * sizeof(float));

  memcpy(Du1, par10+Ndim2, Ndim2 * sizeof(float));
  memcpy(Du1+Ndim2, par11+Ndim2, Ndim2 * sizeof(float));

  Vc::float_v Su0Vec, Du0Vec, Su1Vec, Du1Vec, parUVec;
  float parU[12];
  std::vector<Vc::float_v> parUVector(4);
  for(size_t i = 0; i < Ndim4; i += Vc::float_v::size())
  {
    Su0Vec = Vc::float_v(Su0+i, Vc::Aligned); Du0Vec = Vc::float_v(Du0+i, Vc::Aligned); Su1Vec = Vc::float_v(Su1+i, Vc::Aligned); Du1Vec = Vc::float_v(Du1+i, Vc::Aligned);
    gridU.interpolate(knotU, Su0Vec, Du0Vec, Su1Vec, Du1Vec, u, parUVec);
    parUVec.store(parU+i, Vc::Aligned);
  }
  for(size_t i = 0; i < 4; i++)
  {
    parUVector[i] = Vc::float_v(parU + i * 3);
  }
  Vc::float_v SuVec;
  gridV.interpolate(knotV, parUVector[0], parUVector[1], parUVector[2], parUVector[3], v, SuVec);
  Suv[0] = SuVec[0];
  Suv[1] = SuVec[1];
  Suv[2] = SuVec[2];
}

GPUdi() void Spline2D::interpolateVecVertical(int Ndim, GPUgeneric() const float* parameters, float u, float v, GPUgeneric() float Suv[/*Ndim*/]) const
{
  // Get interpolated value for f(u,v) using parameters[getNumberOfParameters()]
  const Spline1D& gridU = getGridU();
  const Spline1D& gridV = getGridV();
  int nu = gridU.getNumberOfKnots();
  int iu = gridU.getKnotIndex(u);
  int iv = gridV.getKnotIndex(v);

  const Spline1D::Knot& knotU = gridU.getKnot(iu);
  const Spline1D::Knot& knotV = gridV.getKnot(iv);

  const int Ndim2 = Ndim * 2; 
  const int Ndim4 = Ndim * 4; 
  
  // X:=Sx, Y:=Sy, Z:=Sz
  const float* par00 = parameters + (nu * iv + iu) * Ndim4; // values { {X,Y,Z}, {X,Y,Z}'v, {X,Y,Z}'u, {X,Y,Z}''vu } at {u0, v0}
  const float* par10 = par00 + Ndim4;                       // values { ... } at {u1, v0}
  const float* par01 = par00 + Ndim4 * nu;                  // values { ... } at {u0, v1}
  const float* par11 = par01 + Ndim4;                       // values { ... } at {u1, v1}

  Vc::float_v Su00(par00), Su01(par00+3), Su02(par01), Su03(par01+3); 
  Vc::float_v Du00(par00+6), Du01(par00+9), Du02(par01+6), Du03(par01+9); 
  Vc::float_v Su10(par10), Su11(par10+3), Su12(par11), Su13(par11+3); 
  Vc::float_v Du10(par10+6), Du11(par10+9), Du12(par11+6), Du13(par11+9); 
  Vc::float_v parU0, parU1, parU2, parU3, parV; 
  gridU.interpolate(knotU, Su00, Du00, Su10, Du10, u, parU0);
  gridU.interpolate(knotU, Su01, Du01, Su11, Du11, u, parU1);
  gridU.interpolate(knotU, Su02, Du02, Su12, Du12, u, parU2);
  gridU.interpolate(knotU, Su03, Du03, Su13, Du13, u, parU3);
  gridV.interpolate(knotV, parU0, parU1, parU2, parU3, v, parV);
  Suv[0] = parV[0];
  Suv[1] = parV[1];
  Suv[2] = parV[2];
}

GPUdi() void Spline2D::interpolateVecHorizontal(int Ndim, GPUgeneric() const float* parameters, Vc::float_v u, Vc::float_v v, GPUgeneric() Vc::float_v Suv[/*Ndim*/]) const
{
  const Spline1D& gridU = getGridU();
  const Spline1D& gridV = getGridV();
  int nu = gridU.getNumberOfKnots();
  auto iu = gridU.getKnotIndex(u);
  auto iv = gridV.getKnotIndex(v);
  auto knotU = gridU.getKnot(iu);  
  auto knotV = gridV.getKnot(iv);  
  const int Ndim2 = Ndim * 2; 
  const int Ndim4 = Ndim * 4; 
  
  auto par00 = (nu * iv + iu) * Ndim4;  // values { {X,Y,Z}, {X,Y,Z}'v, {X,Y,Z}'u, {X,Y,Z}''vu } at {u0, v0}
  auto par10 = par00 + Ndim4;           // values { ... } at {u1, v0}
  auto par01 = par00 + Ndim4 * nu;      // values { ... } at {u0, v1}
  auto par11 = par01 + Ndim4;           // values { ... } at {u1, v1}   
  
  Vc::float_v Su0[Ndim4], Du0[Ndim4], Su1[Ndim4], Du1[Ndim4], parU[Ndim4];
  for(size_t i = 0; i < Ndim2; i++)
  {
    for(size_t j = 0; j < Vc::float_v::size(); j++)
    {
      Su0[i][j] = *(parameters + par00[j] + i); 
      Su0[Ndim2 + i][j] = *(parameters + par01[j] + i);

      Du0[i][j] = *(parameters + par00[j] + Ndim2 + i);
      Du0[Ndim2 + i][j] = *(parameters + par01[j] + Ndim2 + i);

      Su1[i][j] = *(parameters + par10[j] + i);
      Su1[Ndim2 + i][j] = *(parameters + par11[j] + i);

      Du1[i][j] = *(parameters + par10[j] + Ndim2 + i);
      Du1[Ndim2 + i][j] = *(parameters + par11[j] + Ndim2 + i);
    }
  }
  gridU.interpolate<Vc::float_v>(Ndim4, knotU, Su0, Du0, Su1, Du1, u, parU);
  gridV.interpolate<Vc::float_v>(Ndim, knotV, parU, parU + 3, parU + 6, parU + 9, v, Suv);
}

GPUdi() void Spline2D::interpolateVecHorizontalV2(int Ndim, GPUgeneric() const float* parameters, Vc::float_v u, Vc::float_v v, GPUgeneric() Vc::float_v Suv[/*Ndim*/]) const
{
  const Spline1D& gridU = getGridU();
  const Spline1D& gridV = getGridV();
  int nu = gridU.getNumberOfKnots();
  auto iu = gridU.getKnotIndex(u);
  auto iv = gridV.getKnotIndex(v);
  auto knotU = gridU.getKnot(iu);  
  auto knotV = gridV.getKnot(iv);  
  const int Ndim2 = Ndim * 2; 
  const int Ndim4 = Ndim * 4; 

  const float *par00[Vc::float_v::size()], *par10[Vc::float_v::size()], *par01[Vc::float_v::size()], *par11[Vc::float_v::size()];  
  float Su0[Vc::float_v::size()][Ndim4], Du0[Vc::float_v::size()][Ndim4], Su1[Vc::float_v::size()][Ndim4], Du1[Vc::float_v::size()][Ndim4];
  Vc::float_v Su0Vec, Du0Vec, Su1Vec, Du1Vec, parU[Ndim4];
  for(size_t i = 0; i < Vc::float_v::size(); i++)
  {
    par00[i] = (nu * iv[i] + iu[i]) * Ndim4 + parameters;
    par10[i] = par00[i] + Ndim4;          
    par01[i] = par00[i] + Ndim4 * nu;      
    par11[i] = par01[i] + Ndim4;           
  }
  for(size_t i = 0; i < Vc::float_v::size(); i++)
  {
    memcpy(Su0[i], par00[i], Ndim2 * sizeof(float));
    memcpy(Su0[i] + Ndim2, par01[i], Ndim2 * sizeof(float));

    memcpy(Du0[i], par00[i] + Ndim2, Ndim2 * sizeof(float));
    memcpy(Du0[i] + Ndim2, par01[i] + Ndim2, Ndim2 * sizeof(float));

    memcpy(Su1[i], par10[i], Ndim2 * sizeof(float));
    memcpy(Su1[i] + Ndim2, par11[i], Ndim2 * sizeof(float));

    memcpy(Du1[i], par10[i] + Ndim2, Ndim2 * sizeof(float));
    memcpy(Du1[i] + Ndim2, par11[i] + Ndim2, Ndim2 * sizeof(float));
  }
  for(size_t i = 0; i < Ndim4; i++)
  {
    for(size_t j = 0; j < Vc::float_v::size(); j++)
    {
      Su0Vec[j] = Su0[j][i];
      Du0Vec[j] = Du0[j][i];
      Su1Vec[j] = Su1[j][i];
      Du1Vec[j] = Du1[j][i];
    }
    gridU.interpolate<Vc::float_v>(knotU, Su0Vec, Du0Vec, Su1Vec, Du1Vec, u, parU[i]);
  }
  gridV.interpolate<Vc::float_v>(Ndim, knotV, parU, parU + 3, parU + 6, parU + 9, v, Suv);
}

GPUdi() void Spline2D::interpolateVecHorizontalV3(int Ndim, GPUgeneric() const float* parameters, Vc::float_v u, Vc::float_v v, GPUgeneric() Vc::float_v Suv[/*Ndim*/]) const
{
  const Spline1D& gridU = getGridU();
  const Spline1D& gridV = getGridV();
  int nu = gridU.getNumberOfKnots();
  auto iu = gridU.getKnotIndex(u);
  auto iv = gridV.getKnotIndex(v);
  auto knotU = gridU.getKnot(iu);  
  auto knotV = gridV.getKnot(iv);  
  const int Ndim2 = Ndim * 2; 
  const int Ndim4 = Ndim * 4; 

  const float *par00[Vc::float_v::size()], *par10[Vc::float_v::size()], *par01[Vc::float_v::size()], *par11[Vc::float_v::size()];  
  for(size_t i = 0; i < Vc::float_v::size(); i++)
  {
    par00[i] = parameters + (nu * iv[i] + iu[i]) * Ndim4;
    par10[i] = par00[i] + Ndim4;          
    par01[i] = par00[i] + Ndim4 * nu;      
    par11[i] = par01[i] + Ndim4;           
  }

  Vc::float_v parU[Ndim4];
  for(size_t i = 0; i < Ndim2; i++)
  {
    gridU.interpolate<Vc::float_v>(knotU, Vc::float_v::IndexesFromZero().apply([&](int index){return *(par00[index] + i);}), 
                                          Vc::float_v::IndexesFromZero().apply([&](int index){return *(par00[index] + Ndim2 + i);}),
                                          Vc::float_v::IndexesFromZero().apply([&](int index){return *(par10[index] + i);}), 
                                          Vc::float_v::IndexesFromZero().apply([&](int index){return *(par10[index] + Ndim2 + i);}), 
                                          u, parU[i]);
    gridU.interpolate<Vc::float_v>(knotU, Vc::float_v::IndexesFromZero().apply([&](int index){return *(par01[index] + i);}), 
                                          Vc::float_v::IndexesFromZero().apply([&](int index){return *(par01[index] + Ndim2 + i);}),
                                          Vc::float_v::IndexesFromZero().apply([&](int index){return *(par11[index] + i);}), 
                                          Vc::float_v::IndexesFromZero().apply([&](int index){return *(par11[index] + Ndim2 + i);}), 
                                          u, parU[i + Ndim2]); 
  }
  gridV.interpolate<Vc::float_v>(Ndim, knotV, parU, parU + 3, parU + 6, parU + 9, v, Suv);
}

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
