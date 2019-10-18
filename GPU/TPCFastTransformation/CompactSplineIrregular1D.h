// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  CompactSplineIrregular1D.h
/// \brief Definition of CompactSplineIrregular1D class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#ifndef ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_COMPACTSPLINEIRREGULAR1D_H
#define ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_COMPACTSPLINEIRREGULAR1D_H

#include "GPUCommonDef.h"
#include "FlatObject.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
///
/// The CompactSplineIrregular1D class represents spline interpolation on an one-dimensional irregular (nonunifom) grid.
///
/// The class is a flat C structure. No virtual methods, no ROOT types are used.
///
/// ---
/// The spline interpolates a generic function F:[0,Umax]->R^m.
///
/// The function parameter is called U, the function value is called F (may be multi-dimensional).
/// The interpolation is performed on n knots {U0==0., U1, .., Un-1==Umax}
/// using the function values Fi and the derivatives Di at the knots.
///
/// --- Knots ---
///
/// Umax is an integer number.
/// The knots have integer coordinates on the segment [0,Umax].
/// It is done this way for fast indexing of the segments between knots.
///
/// To interpolate on an interval other than [0,Umax], one should scale the U coordinate and the derivatives Di.
///
/// --- Function values at knots---
///
/// Nothing which depends on F is stored in the class,
/// therefore one can use the same class for interpolation of different input functions on the same knots.
/// The function values F_i and the derivatives D_i = {F'_u}_i have to be provided by the user for each call.
/// The format of the spline input data: { {Fx,Fy,Fz,Dx,Dy,Dz}_0, ..., {Fx,Fy,Fz,Dx,Dy,Dz}_n-1 } for 3-dimensional F
///
/// --- Interpolation ---
/// An interpolation of F values between the knots is performed by 3-th degree polynoms.
/// The polynoms and they 1-st derivatives are continuous at the knots.
/// Depending on the initialization of the derivatives, the second derivative may or may not be continious.
///
/// ---- Initialisation ---
/// The minimal number of knots is 2, the minimal Umax is 1
///
/// Knot U0=0. is always present. It will be added automatically when it is not set by the user.
/// Attention! The number of knots may change during the initialization.
///
/// The user should provide function values Fi and the derivatives Di for all constructed(!) knots.
/// They can be calculated using utilities from the CompactSplineHelper1D class.
///
/// ------------
///
///  Example of creating a spline:
///
///  const int nKnots=3;
///  int knots[nKnots] = {0, 1, 5};
///  CompactSplineIrregular1D spline;
///  spline.construct(nKnots, knots );
///  {// manual
///    float data[2*nKnots] = { 3.5, 0.01, 2.0, -0.01, 3.1, 0.02};
///    spline.getSpline( data, 0.0 ); // == 3.5
///    spline.getSpline( data, 0.2 ); // == some interpolated value
///    spline.getSpline( data, 1.0 ); // == 2.0
///    spline.getSpline( data, 5.0 ); // == 3.1
///  }
///  { // using helper
///    auto F = [&](float u) -> float {
///     return ...; // F(u)
///    };
///    CompactSplineHelper1D helper;
///    std::unique_ptr<float[]> data = helper.constructSpline(spline, F, 0.f, 5.f, 2);
///    spline.getSpline( data.get(), 0.0 ); // == F(0.0)
///    spline.getSpline( data.get(), 0.2 ); // some interpolated value
///    spline.getSpline( data.get(), 1.0 ); // == F(1.0)
///    spline.getSpline( data.get(), 5.0 ); // == F(5.0)
///   }
///
///  --- See also CompactSplineIrregular1D::test();
///
class CompactSplineIrregular1D : public FlatObject
{
 public:
  ///
  /// \brief The struct represents a knot(i) and interval [ knot(i), knot(i+1) ]
  ///
  struct Knot {
    float u;  ///< u coordinate of the knot i (an integer number in float format)
    float Li; ///< inverse length of the [knot_i, knot_{i+1}] segment ( 1./ a (small) integer number)
  };

  /// _____________  Version control __________________________

  /// Version number
  GPUd() static constexpr int getVersion() { return 1; }

  /// Size of the data array in elements, must be multiplied by sizeof(float)
  template <typename T>
  GPUd() size_t getDataSize(int Ndim = 1) const
  {
    return (2 * Ndim * sizeof(T)) * mNumberOfKnots;
  }

  /// Size of the data array in elements, must be multiplied by sizeof(float)

  GPUd() size_t getDataSizeInElements(int Ndim = 1) const { return (2 * Ndim) * mNumberOfKnots; }

  /// _____________  Constructors / destructors __________________________

  /// Default constructor. Creates an empty uninitialised object
  CompactSplineIrregular1D();

  /// Copy constructor: disabled to avoid ambiguity. Use cloneFromObject instead
  CompactSplineIrregular1D(const CompactSplineIrregular1D&) CON_DELETE;

  /// Assignment operator: disabled to avoid ambiguity. Use cloneFromObject instead
  CompactSplineIrregular1D& operator=(const CompactSplineIrregular1D&) CON_DELETE;

  /// Destructor
  ~CompactSplineIrregular1D() CON_DEFAULT;

  /// _____________  FlatObject functionality, see FlatObject class for description  ____________

  /// Memory alignment

  using FlatObject::getBufferAlignmentBytes;
  using FlatObject::getClassAlignmentBytes;

  /// Construction interface

#if !defined(GPUCA_GPUCODE)
  void cloneFromObject(const CompactSplineIrregular1D& obj, char* newFlatBufferPtr);
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

  /// Get minimal required alignment for the spline data
  static constexpr size_t getDataAlignmentBytes() { return 2 * sizeof(float); }

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

  /// Constructor for a regular spline. Knots will be placed at positions i/(numberOfKnots-1)
  void constructRegular(int numberOfKnots);
#endif

  /// _______________  Main functionality   ________________________

  /// Get interpolated value for {F(u): float -> T^Ndim} at the segment [knotL, next knotR] with function values Fl, Fr and slopes Dl, Dr
  template <typename T, int Ndim = 1>
  GPUd() static void getSpline(const CompactSplineIrregular1D::Knot& knotL,
                               GPUgeneric() const T Fl[], GPUgeneric() const T Dl[],
                               GPUgeneric() const T Fr[], GPUgeneric() const T Dr[],
                               float u, GPUgeneric() T Fu[]);

  /// Get interpolated value for F(u) using spline data with a border check
  template <typename T, int Ndim = 1>
  GPUd() void getSpline(GPUgeneric() const T data[], float u, GPUgeneric() T Fu[]) const;

  /// Get interpolated value for F(u) using spline data with no border check
  template <typename T, int Ndim = 1>
  GPUd() void getSplineNonsafe(GPUgeneric() const T data[], float u, GPUgeneric() T Fu[]) const;

  /// Simple interface for 1D spline
  template <typename T>
  GPUd() static T getSpline(const CompactSplineIrregular1D::Knot& knotL,
                            const T& Fl, const T& Dl, const T& Fr, const T& Dr, float u);

  /// Simple interface for 1D spline
  template <typename T>
  GPUd() T getSpline(GPUgeneric() const T data[], float u) const;

  /// _______________  Getters   ________________________

  /// Get number of knots
  GPUd() int getNumberOfKnots() const { return mNumberOfKnots; }

  /// Get index of associated knot for a given U coordinate. With a border check.
  GPUd() int getKnotIndex(float u) const;

  /// Get index of associated knot for a given U coordinate. No border check.
  GPUd() int getKnotIndexNonsafe(float u) const;

  /// Get i-th knot, no border check performed!
  GPUd() const CompactSplineIrregular1D::Knot& getKnot(int i) const { return getKnots()[i]; }

  /// Get the array of knots
  GPUd() const CompactSplineIrregular1D::Knot* getKnots() const { return reinterpret_cast<const CompactSplineIrregular1D::Knot*>(mFlatBufferPtr); }

  /// technical stuff

  /// Get a map  (U axis bin index) -> (corresponding knot index)
  GPUd() const int* getBin2KnotMap() const { return mBin2KnotMap; }

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
  CompactSplineIrregular1D::Knot* getKnotsNonConst() { return reinterpret_cast<CompactSplineIrregular1D::Knot*>(mFlatBufferPtr); }

  /// Non-const accessor to bins->knots map
  int* getBin2KnotMapNonConst() { return mBin2KnotMap; }

  ///
  /// ====  Data members   ====
  ///

  int mNumberOfKnots; ///< n knots on the grid
  int mUmax;          ///< number of axis bins
  int* mBin2KnotMap;  ///< pointer to (axis bin) -> (knot) map inside the mFlatBufferPtr array
};

/// ====================================================
///       Inline implementations of some methods
/// ====================================================

GPUdi() int CompactSplineIrregular1D::getKnotIndexNonsafe(float u) const
{
  /// Get i: u is in [knot_i, knot_{i+1}) interval
  /// no border check! u must be in [0,mUmax]
  return getBin2KnotMap()[(int)u];
}

GPUdi() int CompactSplineIrregular1D::getKnotIndex(float u) const
{
  /// Get i: u is in [knot_i, knot_{i+1}) interval
  /// when u is otside of [0, mUmax], return the edge intervals
  int ibin = (int)u;
  if (ibin < 0)
    ibin = 0;
  if (ibin > mUmax)
    ibin = mUmax;
  return getBin2KnotMap()[ibin];
}

template <typename T, int Ndim = 1>
GPUdi() void CompactSplineIrregular1D::getSpline(const CompactSplineIrregular1D::Knot& knotL,
                                                 GPUgeneric() const T Fl[], GPUgeneric() const T Dl[],
                                                 GPUgeneric() const T Fr[], GPUgeneric() const T Dr[],
                                                 float u, GPUgeneric() T Fu[])
{
  /// A static method.
  /// Gives interpolated value of N-dimensional F(u) at u
  /// input: Fl,Dl,Fr,Dr[Ndim] - N-dim function values and slopes at knots {knotL,knotR}
  /// output: Fu[Ndim] - N-dim interpolated value for F(u)

  T uu = T(u - knotL.u);
  T li = T(knotL.Li);
  T x = uu * li; // scaled u
  if
    constexpr(Ndim == 1)
    { // help an optimizer to not loop over dimensions when N dimensions = 1
      T df = ((*Fr) - (*Fl)) * li;
      T a = (*Dl) + (*Dr) - df - df;
      T b = df - (*Dl) - a;
      *Fu = ((a * x + b) * x + (*Dl)) * uu + (*Fl);
    }
  else {
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

template <typename T, int Ndim = 1>
GPUdi() void CompactSplineIrregular1D::getSpline(GPUgeneric() const T data[], float u, GPUgeneric() T Fu[]) const
{
  /// Get interpolated value for F(u) using data array data[Ndim*2*getNumberOfKnots()].
  /// data = { {Fx,Fy,Fz,Dx,Dy,Dz}_0, ... ,{Fx,Fy,Fz,Dx,Dy,Dz}_{n-1} } for f:u->{x,y,z} case
  /// Safe calculation of the knot index
  int iknot = getKnotIndex(u);
  GPUgeneric() const T* d = data + (2 * Ndim) * iknot;
  getSpline<T, Ndim>(getKnot(iknot), &(d[0]), &(d[Ndim]), &(d[2 * Ndim]), &(d[3 * Ndim]), u, Fu);
}

template <typename T, int Ndim = 1>
GPUdi() void CompactSplineIrregular1D::getSplineNonsafe(GPUgeneric() const T data[], float u, GPUgeneric() T Fu[]) const
{
  /// Get interpolated value for f(u) using data array data[Ndim*2*getNumberOfKnots()].
  /// data = { {Fx,Fy,Fz,Dx,Dy,Dz}_0, ... ,{Fx,Fy,Fz,Dx,Dy,Dz}_{n-1} } for f:u->{x,y,z} case
  /// Non-safe calculation of the knot index.
  int iknot = getKnotIndexNonsafe(u);
  GPUgeneric() const T* d = data + (2 * Ndim) * iknot;
  return getSpline<T, Ndim>(getKnot(iknot), d[0], d[Ndim], d[2 * Ndim], d[3 * Ndim], u, Fu);
}

template <typename T>
GPUdi() T CompactSplineIrregular1D::getSpline(const CompactSplineIrregular1D::Knot& knotL,
                                              const T& Fl, const T& Dl,
                                              const T& Fr, const T& Dr,
                                              float u)
{
  /// Simple interface for 1D spline
  T Fu;
  getSpline(knotL, &Fl, &Dl, &Fr, &Dr, u, &Fu);
  return Fu;
}

template <typename T>
GPUdi() T CompactSplineIrregular1D::getSpline(GPUgeneric() const T data[], float u) const
{
  /// Simple interface for 1D spline
  T Fu;
  getSpline<T, 1>(data, u, &Fu);
  return Fu;
}

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif