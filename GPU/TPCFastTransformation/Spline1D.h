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
/// The interpolation is performed on n knots {U0==0., U1, .., Un-1==Umax}
/// using the function values Fi and the derivatives Di at the knots.
///
/// --- Knots ---
///
/// Umax is an integer number.
/// The knots belong to [0, Umax] and have integer coordinates.
/// It is implemented this way for fast matching of an U coordinate to the corresponding segment between the knots.
///
/// To interpolate a function on an interval other than [0,Umax], one should scale the U coordinate and the derivatives Di.
///
/// --- Function values at knots---
///
/// Nothing which depends on F is stored in the class,
/// therefore one can use the same class for interpolation of different input functions on the same knots.
/// The function values F_i and the derivatives D_i == {F'_u}_i have to be provided by the user for each call.
/// The format of the spline data: { {Fx,Fy,Fz,Dx,Dy,Dz}_0, ..., {Fx,Fy,Fz,Dx,Dy,Dz}_n-1 } for 3-dimensional F
///
/// --- Interpolation ---
///
/// The interpolation of F values between the knots is performed by the 3-th degree polynoms.
/// The polynoms and they 1-st derivatives are continuous at the knots.
/// Depending on the initialization of the derivatives, the second derivative may or may not be continious.
///
/// ---- Initialisation of the spline knots ---
///
/// The minimal number of knots is 2, the minimal Umax is 1
///
/// The knot U0=0.0 must always be present. It will be added automatically when it is not set by the user.
/// Duplicated knots are not allowed and will be removed.
/// Therefore: Attention! The number of knots may change during the initialization.
///
/// ---- Initialisation of the spline data for a given F ---
///
/// The user should create an array of the function values Fi and the derivatives Di at all (initialised!) knots.
/// This data array can be created using utilities from the SplineHelper1D class.
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
///    float data[2*nKnots] = { F(0.), 0.01, F(2.), -0.01, F(3.), 0.02};
///    spline.getSpline<1>( data, 0.0 ); // == F(0.)
///    spline.getSpline<1>( data, 0.2 ); // == some interpolated value
///    spline.getSpline<1>( data, 1.0 ); // == F(2.)
///    spline.getSpline<1>( data, 5.0 ); // == F(3.)
///  }
///  { // using the helper
///    SplineHelper1D helper;
///    std::unique_ptr<float[]> data = helper.constructData(spline, F, 0.f, 1.f, 2);
///    spline.getSpline<1>( data.get(), 0.0 ); // == F(0.)
///    spline.getSpline<1>( data.get(), 0.2 ); // some interpolated value
///    spline.getSpline<1>( data.get(), 1.0 ); // == F(0.2)
///    spline.getSpline<1>( data.get(), 5.0 ); // == F(1.)
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

  /// Default constructor. Creates an empty uninitialised object
  Spline1D();

  /// Copy constructor: disabled to avoid ambiguity. Use cloneFromObject instead
  Spline1D(const Spline1D&) CON_DELETE;

  /// Assignment operator: disabled to avoid ambiguity. Use cloneFromObject instead
  Spline1D& operator=(const Spline1D&) CON_DELETE;

  /// Destructor
  ~Spline1D() CON_DEFAULT;

  /// _____________  FlatObject functionality, see FlatObject class for description  ____________

  /// Memory alignment

  using FlatObject::getBufferAlignmentBytes;
  using FlatObject::getClassAlignmentBytes;

  /// Construction interface

#if !defined(GPUCA_GPUCODE)
  void cloneFromObject(const Spline1D& obj, char* newFlatBufferPtr);
#endif
  void destroy();

  /// Making the data buffer external

  using FlatObject::releaseInternalBuffer;

#if !defined(GPUCA_GPUCODE)
  void moveBufferTo(char* newBufferPtr);
#endif

  /// Moving the class with its external buffer to another location

  void setActualBufferAddress(char* actualFlatBufferPtr);
  void setFutureBufferAddress(char* futureFlatBufferPtr);

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

  /// Get interpolated value for {F(u): float -> T^Ndim} at the segment [knotL, next knotR] with function values Fl, Fr and slopes Dl, Dr
  template <int Ndim, typename T>
  GPUd() static void getSpline(const Spline1D::Knot& knotL,
                               GPUgeneric() const T Fl[Ndim], GPUgeneric() const T Dl[Ndim],
                               GPUgeneric() const T Fr[Ndim], GPUgeneric() const T Dr[Ndim],
                               float u, GPUgeneric() T Fu[Ndim]);

  /// Get interpolated value for F(u) using spline data with a border check
  template <int Ndim, typename T>
  GPUd() void getSpline(GPUgeneric() const T data[], float u, GPUgeneric() T Fu[Ndim]) const;

  /// Get interpolated value for F(u) using spline data with no border check
  template <int Ndim, typename T>
  GPUd() void getSplineNonSafe(GPUgeneric() const T data[], float u, GPUgeneric() T Fu[Ndim]) const;

  /// Simple interface for 1D spline
  template <typename T>
  GPUd() static T getSpline(GPUgeneric() const Spline1D::Knot& knotL,
                            GPUgeneric() const T& Fl, GPUgeneric() const T& Dl, GPUgeneric() const T& Fr, GPUgeneric() const T& Dr, float u);

  /// Simple interface for 1D spline
  template <typename T>
  GPUd() T getSpline(GPUgeneric() const T data[], float u) const;

  /// _______________  Getters   ________________________

  /// Get minimal required alignment for the spline data
  template <int Ndim, typename T>
  static constexpr size_t getDataAlignmentBytes()
  {
    return std::min<2 * sizeof(T) * Ndim, 16>;
  }

  /// Size of the data array in bytes
  template <int Ndim, typename T>
  GPUd() size_t getDataSize() const
  {
    return (2 * Ndim * sizeof(T)) * mNumberOfKnots;
  }

  /// Size of the data array in elements
  template <int Ndim>
  GPUd() size_t getDataSizeInElements() const
  {
    return (2 * Ndim) * mNumberOfKnots;
  }

  /// Get number of knots
  GPUd() int getNumberOfKnots() const { return mNumberOfKnots; }

  /// Get index of associated knot for a given U coordinate. With a border check.
  GPUd() int getKnotIndex(float u) const;

  /// Get index of associated knot for a given U coordinate. No border check.
  GPUd() int getKnotIndexNonSafe(float u) const;

  /// Get i-th knot, no border check performed!
  GPUd() const Spline1D::Knot& getKnot(int i) const { return getKnots()[i]; }

  /// Get the array of knots
  GPUd() const Spline1D::Knot* getKnots() const { return reinterpret_cast<const Spline1D::Knot*>(mFlatBufferPtr); }

  /// technical stuff

  /// Get a map (integer U -> corresponding knot index)
  GPUd() const int* getUtoKnotMap() const { return mUtoKnotMap; }

  /// Get number of axis bins
  int getUmax() const { return mUmax; }

  /// Print method
  void print() const;

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE) // code invisible on GPU and in the standalone compilation
  /// Test the class functionality
  static int test(bool draw = 0);
#endif

 private:
  /// Non-const accessor to knots array
  Spline1D::Knot* getKnotsNonConst() { return reinterpret_cast<Spline1D::Knot*>(mFlatBufferPtr); }

  /// Non-const accessor to U->knots map
  int* getUtoKnotMapNonConst() { return mUtoKnotMap; }

  ///
  /// ====  Data members   ====
  ///

  int mNumberOfKnots; ///< n knots on the grid
  int mUmax;          ///< number of axis bins
  int* mUtoKnotMap;   ///< pointer to (integer U -> knot index) map inside the mFlatBufferPtr array
};

/// ====================================================
///       Inline implementations of some methods
/// ====================================================

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

template <int Ndim, typename T>
GPUdi() void Spline1D::getSpline(const Spline1D::Knot& knotL,
                                        GPUgeneric() const T Fl[Ndim], GPUgeneric() const T Dl[Ndim],
                                        GPUgeneric() const T Fr[Ndim], GPUgeneric() const T Dr[Ndim],
                                        float u, GPUgeneric() T Fu[Ndim])
{
  /// A static method.
  /// Gives interpolated value of N-dimensional F(u) at u
  /// input: Fl,Dl,Fr,Dr[Ndim] - N-dim function values and slopes at knots {knotL,knotR}
  /// output: Fu[Ndim] - N-dim interpolated value for F(u)

  T uu = T(u - knotL.u);
  T li = T(knotL.Li);
  T x = uu * li;             // scaled u
  if constexpr (Ndim == 1) { // help the optimizer to not loop over dimensions when Ndim == 1
    T df = ((*Fr) - (*Fl)) * li;
    T a = (*Dl) + (*Dr) - df - df;
    T b = df - (*Dl) - a;
    *Fu = ((a * x + b) * x + (*Dl)) * uu + (*Fl);
  } else {
    for (int i = 0; i < Ndim; ++i) {
      T df = (Fr[i] - Fl[i]) * li;
      T a = Dl[i] + Dr[i] - df - df;
      T b = df - Dl[i] - a;
      Fu[i] = ((a * x + b) * x + Dl[i]) * uu + Fl[i];
    }
  }

  /* another way to calculate f(u):
  T uu = T(u - knotL.u);
  T x = uu * T(knotL.Li); // scaled u
  T xm1 = x-1;
  T x2 = x * x;
  float cFr = x2*(3-2*x);
  float cFl = 1-cFr;
  float cDl = x*xm1*xm1*knotL.L;
  float cDr = x2*xm1*knotL.L;
  return cFl*Fl + cFr*Fr + cDl*Dl + cDr*Dr;
  */
} // namespace gpu

template <int Ndim, typename T>
GPUdi() void Spline1D::getSpline(GPUgeneric() const T data[], float u, GPUgeneric() T Fu[Ndim]) const
{
  /// Get interpolated value for F(u) using data array data[Ndim*2*getNumberOfKnots()].
  /// data = { {Fx,Fy,Fz,Dx,Dy,Dz}_0, ... ,{Fx,Fy,Fz,Dx,Dy,Dz}_{n-1} } for f:u->{x,y,z} case
  /// Safe calculation of the knot index
  int iknot = getKnotIndex(u);
  const T* d = data + (2 * Ndim) * iknot;
  getSpline<Ndim, T>(getKnot(iknot), &(d[0]), &(d[Ndim]), &(d[2 * Ndim]), &(d[3 * Ndim]), u, Fu);
}

template <int Ndim, typename T>
GPUdi() void Spline1D::getSplineNonSafe(GPUgeneric() const T data[], float u, GPUgeneric() T Fu[]) const
{
  /// Get interpolated value for f(u) using data array data[Ndim*2*getNumberOfKnots()].
  /// data = { {Fx,Fy,Fz,Dx,Dy,Dz}_0, ... ,{Fx,Fy,Fz,Dx,Dy,Dz}_{n-1} } for f:u->{x,y,z} case
  /// Non-safe calculation of the knot index.
  int iknot = getKnotIndexNonSafe(u);
  GPUgeneric() const T* d = data + (2 * Ndim) * iknot;
  return getSpline<Ndim, T>(getKnot(iknot), d[0], d[Ndim], d[2 * Ndim], d[3 * Ndim], u, Fu);
}

template <typename T>
GPUdi() T Spline1D::getSpline(GPUgeneric() const Spline1D::Knot& knotL,
                                     GPUgeneric() const T& Fl, GPUgeneric() const T& Dl,
                                     GPUgeneric() const T& Fr, GPUgeneric() const T& Dr,
                                     float u)
{
  /// Simple interface for 1D spline
  T Fu;
  getSpline<1>(knotL, &Fl, &Dl, &Fr, &Dr, u, &Fu);
  return Fu;
}

template <typename T>
GPUdi() T Spline1D::getSpline(GPUgeneric() const T data[], float u) const
{
  /// Simple interface for 1D spline
  T Fu;
  getSpline<1>(data, u, &Fu);
  return Fu;
}

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif