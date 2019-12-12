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
  GPUd() void interpolate(int Ndim, GPUgeneric() const T* parameters, float u, float v, GPUgeneric() T Fuv[/*Ndim*/]) const;

  /// Same as interpolate, but using vectorized calculation.
  /// \param parameters should be at least 128-bit aligned
  template <typename T>
  GPUd() void interpolateVec(int Ndim, GPUgeneric() const T* parameters, float u, float v, GPUgeneric() T Fuv[/*Ndim*/]) const;

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
    return sizeof(T) * getNumberOfParameters(Ndim);
  }

  /// Number of parameters
  GPUd() size_t getNumberOfParameters(int Ndim) const
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
GPUdi() void Spline2D::interpolate(int Ndim, GPUgeneric() const T* parameters, float u, float v, GPUgeneric() T Fuv[/*Ndim*/]) const
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

  const T* parameters00 = parameters + (nu * iv + iu) * Ndim4; // values { {X,Y,Z}, {X,Y,Z}'v, {X,Y,Z}'u, {X,Y,Z}''vu } at {u0, v0}
  const T* parameters10 = parameters00 + Ndim4;                // values { ... } at {u1, v0}
  const T* parameters01 = parameters00 + Ndim4 * nu;           // values { ... } at {u0, v1}
  const T* parameters11 = parameters01 + Ndim4;                // values { ... } at {u1, v1}

  T Fu0[Ndim4]; // values { {X,Y,Z,X'v,Y'v,Z'v}(v0), {X,Y,Z,X'v,Y'v,Z'v}(v1) }, at u0
  T Du0[Ndim4]; // derivatives {}'_u  at u0
  for (int i = 0; i < Ndim2; i++)
    Fu0[i] = parameters00[i];
  for (int i = 0; i < Ndim2; i++)
    Fu0[Ndim2 + i] = parameters01[i];
  for (int i = 0; i < Ndim2; i++)
    Du0[i] = parameters00[Ndim2 + i];
  for (int i = 0; i < Ndim2; i++)
    Du0[Ndim2 + i] = parameters01[Ndim2 + i];

  T Fu1[Ndim4]; // values { {X,Y,Z,X'v,Y'v,Z'v}(v0), {X,Y,Z,X'v,Y'v,Z'v}(v1) }, at u1
  T Du1[Ndim4]; // derivatives {}'_u  at u1
  for (int i = 0; i < Ndim2; i++)
    Fu1[i] = parameters10[i];
  for (int i = 0; i < Ndim2; i++)
    Fu1[Ndim2 + i] = parameters11[i];
  for (int i = 0; i < Ndim2; i++)
    Du0[i] = parameters10[Ndim2 + i];
  for (int i = 0; i < Ndim2; i++)
    Du0[Ndim2 + i] = parameters11[Ndim2 + i];

  T parametersU[Ndim4]; // interpolated values { {X,Y,Z,X'v,Y'v,Z'v}(v0), {X,Y,Z,X'v,Y'v,Z'v}(v1) } at u
  gridU.interpolate<T>(Ndim4, knotU, Fu0, Du0, Fu1, Du1, u, parametersU);

  T* Fv0 = parametersU + 0;
  T* Dv0 = parametersU + Ndim;
  T* Fv1 = parametersU + Ndim2;
  T* Dv1 = parametersU + Ndim2 + Ndim;

  gridV.interpolate<T>(Ndim, knotV, Fv0, Dv0, Fv1, Dv1, v, Fuv);
}

template <typename T>
GPUdi() void Spline2D::interpolateVec(int Ndim, GPUgeneric() const T* parameters, float u, float v, GPUgeneric() T Fuv[/*Ndim*/]) const
{
  // Same as interpolate, but using vectorized calculation.
  // \param parameters should be at least 128-bit aligned

  /// TODO: vectorize
  interpolate<T>(Ndim, parameters, u, v, Fuv);
}

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
