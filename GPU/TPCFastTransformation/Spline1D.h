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

namespace GPUCA_NAMESPACE
{
namespace gpu
{
///
/// The Spline1D class represents spline interpolation on an one-dimensional irregular (nonunifom) grid.
///
/// The class is a flat C structure. No virtual methods, no ROOT types are used.
///
/// ---
/// The spline interpolates a function F:[0,Umax]->R^m.
///
/// The function parameter is called U, the function value is called F. The F may be multi-dimensional.
/// The spline value is called S(u).
/// The interpolation is performed on n knots {U0==0., U1, .., Un-1==Umax}
/// using the spline values Si and the spline derivatives Di==S'i at the knots.
///
/// --- Knots ---
///
/// Umax is an integer number.
/// The knots belong to [0, Umax] and have integer coordinates.
/// It is implemented this way for fast matching of any U coordinate to the corresponding segment between the knots.
///
/// To interpolate a function on an interval other than [0,Umax], one should scale the U coordinate and the derivatives Di.
/// It is implemented this way in order to avoid a double scaling, since during the interpolation 
/// the U coordinate needs to be scalled to [0,Umax] anyhow.
///
/// --- Function values at knots---
///
/// Nothing which depends on F is stored in the class,
/// therefore one can use the same spline object for interpolation of different input functions.
/// The spline parameters S_i and D_i must be provided by the user for each call.
/// The format of the spline parameters: { {Sx,Sy,Sz,Dx,Dy,Dz}_0, ..., {Sx,Sy,Sz,Dx,Dy,Dz}_n-1 } for a 3-dimensional F
///
/// --- Interpolation ---
///
/// The interpolation of F values between the knots is performed by the 3-th degree polynoms.
/// The polynoms and they 1-st derivatives are continuous at the knots.
/// Depending on the initialization of the derivatives, the second derivative may or may not be continious.
///
/// ---- Construction of a spline ---
///
/// The minimal number of knots is 2, the minimal Umax is 1
///
/// The knot U0=0.0 must always be present. It will be added automatically when it is not set by the user.
/// Duplicated knots are not allowed and will be removed.
/// Therefore: Attention! The number of knots may change during the construction.
///
/// ---- Construction of spline parameters for a given F ---
///
/// The user should create an array of the spline values Si and its derivatives Di at all the knots.
/// This parameter array can be created with the help of SplineHelper1D class.
///
/// ---- Flat Object implementation ----
///
/// The class inherits from the FlatObjects. Copying can only be done using the FlatObject interface.
///
/// ------------
///
///  Example of creating a spline:
///
///  const int nKnots=3;
///  int knots[nKnots] = {0, 1, 5};
///  auto F = [&](float u) -> float {
///   return ...; // F(u)
///  };
///  Spline1D spline;
///  spline.construct(nKnots, knots );
///  {// manual construction
///    float parameters[2*nKnots] = { F(0.), 0.01, F(2.), -0.01, F(3.), 0.02};
///    spline.interpolate1D(parameters, 0.0 ); // == F(0.)
///    spline.interpolate1D(parameters, 0.2 ); // == some interpolated value
///    spline.interpolate1D(parameters, 1.0 ); // == F(2.)
///    spline.interpolate1D(parameters, 5.0 ); // == F(3.)
///  }
///  { // using the helper
///    SplineHelper1D helper;
///    helper.setSpline(spline, 2);
///    std::unique_ptr<float[]> parameters = helper.constructParameters(F, 0.f, 1.f);
///    spline.interpolate1D(parameters.get(), 0.0 ); // == F(0.)
///    spline.interpolate1D(parameters.get(), 0.2 ); // some interpolated value
///    spline.interpolate1D(parameters.get(), 1.0 ); // == F(0.2)
///    spline.interpolate1D(parameters.get(), 5.0 ); // == F(1.)
///   }
///
///  --- See also Spline1D::test();
///
class Spline1D : public FlatObject
{
 public:
  ///
  /// \brief The struct Knot represents a knot(i) and the interval [ knot(i), knot(i+1) ]
  ///
  struct Knot {
    float u;  ///< u coordinate of the knot i (an integer number in float format)
    float Li; ///< inverse length of the [knot_i, knot_{i+1}] segment ( == 1.f/ a (small) integer number)
  };

  /// _____________  Version control __________________________

  /// Version number
  GPUd() static constexpr int getVersion() { return 1; }

    /// _____________  Constructors / destructors __________________________

#if !defined(GPUCA_GPUCODE)
  /// Default constructor. Creates a spline with 2 knots.
  Spline1D() : FlatObject() { constructRegular(2); }

  /// Constructor for an irregular spline.
  Spline1D(int numberOfKnots, const int knots[]) : FlatObject() { construct(numberOfKnots, knots); }

  /// Constructor for a regular spline.
  Spline1D(int numberOfKnots) : FlatObject() { constructRegular(numberOfKnots); }
#else
  /// Disable the constructor for the GPU implementation
  Spline1D() CON_DELETE;
#endif

  /// Copy constructor: disabled to avoid ambiguity. Use cloneFromObject instead
  Spline1D(const Spline1D&) CON_DELETE;

  /// Assignment operator: disabled to avoid ambiguity. Use cloneFromObject instead
  Spline1D& operator=(const Spline1D&) CON_DELETE;

  /// Destructor
  ~Spline1D() CON_DEFAULT;

  /// _______________  Construction interface  ________________________

#if !defined(GPUCA_GPUCODE)
  /// Constructor
  ///
  /// Number of created knots may differ from the input values:
  /// - Edge knots {0} and {Umax} will be added if they are not present.
  /// - Duplicated knots, knots with a negative coordinate will be deleted
  /// - At least 2 knots will be created
  ///
  /// \param numberOfKnots     Number of knots in knots[] array
  /// \param knots             Array of knot positions (integer values)
  ///
  void construct(int numberOfKnots, const int knots[]);

  /// Constructor for a regular spline. Knots will be placed at the positions i/(numberOfKnots-1)
  void constructRegular(int numberOfKnots);
#endif

  /// _______________  Main functionality   ________________________

  /// Get interpolated value for {S(u): float -> T^Ndim} at the segment [knotL, next knotR] with the spline values Sl, Sr and the slopes Dl, Dr
  template <typename T>
  GPUd() static void interpolate(int Ndim, const Spline1D::Knot& knotL,
                                 GPUgeneric() const T Sl[/*Ndim*/], GPUgeneric() const T Dl[/*Ndim*/],
                                 GPUgeneric() const T Sr[/*Ndim*/], GPUgeneric() const T Dr[/*Ndim*/],
                                 float u, GPUgeneric() T Su[/*Ndim*/]);

  /// Get interpolated value for F(u) using spline parameters with a border check
  template <typename T>
  GPUd() void interpolate(int Ndim, GPUgeneric() const T parameters[], float u, GPUgeneric() T Su[/*Ndim*/]) const;

  /// Get interpolated value for F(u) using spline parameters with no border check
  template <typename T>
  GPUd() void interpolateNonSafe(int Ndim, GPUgeneric() const T parameters[], float u, GPUgeneric() T Su[/*Ndim*/]) const;

  /// Simplified interface for 1D spline
  template <typename T>
  GPUd() static T interpolate1D(GPUgeneric() const Spline1D::Knot& knotL,
                                GPUgeneric() const T& Sl, GPUgeneric() const T& Dl, GPUgeneric() const T& Sr, GPUgeneric() const T& Dr, float u);

  /// Simplified interface for 1D spline
  template <typename T>
  GPUd() T interpolate1D(GPUgeneric() const T parameters[], float u) const;

  /// _______________  Getters   ________________________

  /// Get U coordinate of the last knot
  int getUmax() const { return mUmax; }

  /// Get minimal required alignment for the spline parameters
  template <typename T>
  static constexpr size_t getParameterAlignmentBytes(int Ndim)
  {
    return std::min<2 * sizeof(T) * Ndim, 16>;
  }

  /// Size of the parameter array in bytes
  template <typename T>
  GPUd() size_t getSizeOfParameters(int Ndim) const
  {
    return sizeof(T) * (size_t) getNumberOfParameters(Ndim);
  }

  /// Number of parameters
  GPUd() int getNumberOfParameters(int Ndim) const
  {
    return (2 * Ndim) * mNumberOfKnots;
  }

  /// Get number of knots
  GPUd() int getNumberOfKnots() const { return mNumberOfKnots; }

  /// Get the array of knots
  GPUd() const Spline1D::Knot* getKnots() const { return reinterpret_cast<const Spline1D::Knot*>(mFlatBufferPtr); }

  /// Get i-th knot
  GPUd() const Spline1D::Knot& getKnot(int i) const { return getKnots()[i < 0 ? 0 : (i >= mNumberOfKnots ? mNumberOfKnots - 1 : i)]; }

  /// Get index of an associated knot for a given U coordinate. Performs a border check.
  GPUd() int getKnotIndex(float u) const;

  /// _______________  Getters with no border check   ________________________

  /// Get index of an associated knot for a given U coordinate. No border check preformed!
  GPUd() int getKnotIndexNonSafe(float u) const;

  /// Get i-th knot. No border check performed!
  GPUd() const Spline1D::Knot& getKnotNonSafe(int i) const { return getKnots()[i]; }

  /// _______________  Technical stuff   ________________________

  /// Get a map (integer U -> corresponding knot index)
  GPUd() const int* getUtoKnotMap() const { return mUtoKnotMap; }

  /// Print method
  void print() const;

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE) // code invisible on GPU and in the standalone compilation
  /// Test the class functionality
  static int test(bool draw = 0);
#endif

  /// _____________  FlatObject functionality, see FlatObject class for description  ____________

  /// Memory alignment

  using FlatObject::getBufferAlignmentBytes;
  using FlatObject::getClassAlignmentBytes;

  /// Construction interface

#if !defined(GPUCA_GPUCODE)
  void cloneFromObject(const Spline1D& obj, char* newFlatBufferPtr);
#endif
  void destroy();

  /// Making the parameter buffer external

  using FlatObject::releaseInternalBuffer;

#if !defined(GPUCA_GPUCODE)
  void moveBufferTo(char* newBufferPtr);
#endif

  /// Moving the class with its external buffer to another location

  void setActualBufferAddress(char* actualFlatBufferPtr);
  void setFutureBufferAddress(char* futureFlatBufferPtr);

 private:
  /// Non-const accessor to knots array
  Spline1D::Knot* getKnotsNonConst() { return reinterpret_cast<Spline1D::Knot*>(mFlatBufferPtr); }

  /// Non-const accessor to U->knots map
  int* getUtoKnotMapNonConst() { return mUtoKnotMap; }

  /// _____________  Data members  ____________

  int mNumberOfKnots; ///< n knots on the grid
  int mUmax;          ///< U of the last knot
  int* mUtoKnotMap;   ///< pointer to (integer U -> knot index) map inside the mFlatBufferPtr array
};

///
/// ========================================================================================================
///       Inline implementations of some methods
/// ========================================================================================================
///

GPUdi() int Spline1D::getKnotIndexNonSafe(float u) const
{
  /// Get i: u is in [knot_i, knot_{i+1}) interval
  /// no border check! u must be in [0,mUmax]
  return getUtoKnotMap()[(int)u];
}

GPUdi() int Spline1D::getKnotIndex(float u) const
{
  /// Get i: u is in [knot_i, knot_{i+1}) interval
  /// when u is otside of [0, mUmax], return the edge intervals
  int iu = (int)u;
  if (iu < 0)
    iu = 0;
  if (iu > mUmax)
    iu = mUmax;
  return getUtoKnotMap()[iu];
}

template <typename T>
GPUdi() void Spline1D::interpolate(int Ndim, const Spline1D::Knot& knotL,
                                   GPUgeneric() const T Sl[/*Ndim*/], GPUgeneric() const T Dl[/*Ndim*/],
                                   GPUgeneric() const T Sr[/*Ndim*/], GPUgeneric() const T Dr[/*Ndim*/],
                                   float u, GPUgeneric() T Su[/*Ndim*/])
{
  /// A static method.
  /// Gives interpolated value of N-dimensional S(u) at u
  /// input: Sl,Dl,Sr,Dr[Ndim] - N-dim function values and slopes at knots {knotL,knotR}
  /// output: Su[Ndim] - N-dim interpolated value for S(u)

  T uu = T(u - knotL.u);
  T li = T(knotL.Li);
  T x = uu * li; // scaled u
  for (int dim = 0; dim < Ndim; ++dim) {
    T df = (Sr[dim] - Sl[dim]) * li;
    T a = Dl[dim] + Dr[dim] - df - df;
    T b = df - Dl[dim] - a;
    Su[dim] = ((a * x + b) * x + Dl[dim]) * uu + Sl[dim];
  }

  /* another way to calculate f(u):
  T uu = T(u - knotL.u);
  T x = uu * T(knotL.Li); // scaled u
  T xm1 = x-1;
  T x2 = x * x;
  float cSr = x2*(3-2*x);
  float cSl = 1-cSr;
  float cDl = x*xm1*xm1*knotL.L;
  float cDr = x2*xm1*knotL.L;
  return cSl*Sl + cSr*Sr + cDl*Dl + cDr*Dr;
  */
}

template <typename T>
GPUdi() void Spline1D::interpolate(int Ndim, GPUgeneric() const T parameters[], float u, GPUgeneric() T Su[/*Ndim*/]) const
{
  /// Get interpolated value for S(u) using parameters[Ndim*2*getNumberOfKnots()].
  /// parameters = { {Sx,Sy,Sz,Dx,Dy,Dz}_0, ... ,{Sx,Sy,Sz,Dx,Dy,Dz}_{n-1} } for f:u->{x,y,z} case
  /// Safe calculation of the knot index
  int iknot = getKnotIndex(u);
  const T* d = parameters + (2 * iknot)*Ndim;
  interpolate<T>(Ndim, getKnotNonSafe(iknot), &(d[0]), &(d[Ndim]), &(d[2 * Ndim]), &(d[3 * Ndim]), u, Su);
}

template <typename T>
GPUdi() void Spline1D::interpolateNonSafe(int Ndim, GPUgeneric() const T parameters[], float u, GPUgeneric() T Su[/*Ndim*/]) const
{
  /// Get interpolated value for f(u) using parameters[Ndim*2*getNumberOfKnots()].
  /// parameters = { {Sx,Sy,Sz,Dx,Dy,Dz}_0, ... ,{Sx,Sy,Sz,Dx,Dy,Dz}_{n-1} } for f:u->{x,y,z} case
  /// Non-safe calculation of the knot index.
  int iknot = getKnotIndexNonSafe(u);
  GPUgeneric() const T* d = parameters + (2 * Ndim) * iknot;
  return interpolate<T>(Ndim, getKnotNonSafe(iknot), d[0], d[Ndim], d[2 * Ndim], d[3 * Ndim], u, Su);
}

template <typename T>
GPUdi() T Spline1D::interpolate1D(GPUgeneric() const Spline1D::Knot& knotL,
                                  GPUgeneric() const T& Sl, GPUgeneric() const T& Dl,
                                  GPUgeneric() const T& Sr, GPUgeneric() const T& Dr,
                                  float u)
{
  /// Simplified interface for 1D spline
  T Su;
  interpolate<T>(1, knotL, &Sl, &Dl, &Sr, &Dr, u, &Su);
  return Su;
}

template <typename T>
GPUdi() T Spline1D::interpolate1D(GPUgeneric() const T parameters[], float u) const
{
  /// Simplified interface for 1D spline
  T Su;
  interpolate<T>(1, parameters, u, &Su);
  return Su;
}

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif