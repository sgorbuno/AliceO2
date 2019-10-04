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
/// The CompactSplineIrregular1D class represents one-dimensional spline interpolation on nonunifom (irregular) grid.
///
/// The class is flat C structure. No virtual methods, no ROOT types are used.
/// It is designed for spline parameterisation of TPC transformation.
///
/// ---
/// The spline interpolates a generic function F:[0,Umax]->R.
///
/// The function parameter is called U, the function value is F.
/// The interpolation is performed on n knots {U0==0., U1, .., Un-1==Umax}
/// using given function values {F0, ..., Fn-1} and derivatives [D0,..,Dn-1] at the knots.
///
/// Umax is an integer number.
/// The knots must have integer coordinates on the segment [0,Umax].
/// It is done this way for a fast indexing of the segments between knots.
///
/// To interpolate on any segment other than [0,Umax], one should scale the U coordinate and the derivatives.
///
/// An interpolation on each segment between two knots is performed by a 3-th degree polynom.
/// The polynoms and they 1-st derivatives are continuous at the knots.
/// Depending on the initialization of the first derivative, the second derivative may or may not be continious.
///
///
/// Nothing which depends on F is stored in the class,
/// therefore one can use the same class for interpolation of different input functions.
/// The function values {F0,..,Fn-1} and the derivatives {D0,..,Dn-1} have to be provided by the user for each call.
///
/// The minimal number of knots is 2, the minimal Umax is 1
///
/// Knot U0=0. is always present. It will be added automatically when it is not set by the user.
/// Attention! The number of knots may change during the initialization.
///
/// The user should provide function values Fi and the derivatives Di for all constructed(!) knots.
/// They can be calculated using the utilities from CompactSplineHelper class.
///
/// ------------
///
///  Example of creating a spline:
///
///  const int nKnots=3;
///  int knots[nKnots] = {0, 1, 5}; // original knot positions
///  CompactSplineIrregular1D spline;
///  spline.construct(nKnots, knots );
///  float data[2*nKnots] = { 3.5, 0.01, 2.0, -0.01, 3.1, 0.02};
///  spline.getSpline( f, 0.0 ); // == 3.5
///  spline.getSpline( f, 0.2 ); // == some interpolated value
///  spline.getSpline( f, 1.0 ); // == 2.0
///  spline.getSpline( f, 5.0 ); // == 3.1
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
  GPUd() size_t getDataSizeInelements() const { return 2 * mNumberOfKnots; }

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

  /// Get interpolated value for F(u) using spline at the segment [knot,next knot] with function values f0, f1 and derivatives d0, d1
  template <typename T>
  GPUd() static T getSpline(const CompactSplineIrregular1D::Knot& knot, T f0, T d0, T f1, T d1, float u);

  /// Get interpolated value for F(u) using data array data[2*getNumberOfKnots()] == {F0,D0, .., Fn-1,Dn-1}
  template <typename T>
  GPUd() GPUgeneric() T getSpline(GPUgeneric() const T data[], float u) const;

  /// Get interpolated value for F(u) using data array data[2*getNumberOfKnots()] == {F0,D0, .., Fn-1,Dn-1}, with a border check
  template <typename T>
  GPUd() GPUgeneric() T getSplineSafe(GPUgeneric() const T data[], float u) const;

  /// Get number of knots
  GPUd() int getNumberOfKnots() const { return mNumberOfKnots; }

  /// Get index of associated knot for a given U coordinate. No border check.
  GPUd() int getKnotIndex(float u) const;

  /// Get index of associated knot for a given U coordinate. With a border check.
  GPUd() int getKnotIndexSafe(float u) const;

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

template <typename T>
GPUdi() T CompactSplineIrregular1D::getSpline(const CompactSplineIrregular1D::Knot& knot0, T f0, T d0, T f1, T d1, float u)
{
  /// A static method.
  /// Get interpolated value for f(u) using spline at knot "knot0" and function & derivative values at knots {knot_0,knot_1}

  T uu = T(u - knot0.u);
  T x = uu * T(knot0.Li); // scaled u

  /* another way to calculate
  T xm1 = x-1;
  T x2 = x * x;
  float cf1 = x2*(3-2*x);
  float cf0 = 1-cf1;
  float cd0 = x*xm1*xm1*knot0.L;
  float cd1 = x2*xm1*knot0.L;
  return cf0*f0 + cf1*f1 + cd0*d0 + cd1*d1;
  */
  f1 = (f1 - f0) * knot0.Li;
  T a = -f1 - f1 + d0 + d1;
  T b = f1 - d0 - a;
  return ((a * x + b) * x + d0) * uu + f0;
}

template <typename T>
GPUdi() GPUgeneric() T CompactSplineIrregular1D::getSplineSafe(GPUgeneric() const T data[], float u) const
{
  /// Get interpolated value for f(u) using data array data[2*getNumberOfKnots()]
  int iknot = getKnotIndexSafe(u);
  const CompactSplineIrregular1D::Knot& knot = getKnot(iknot);
  const T* d = data + (iknot + iknot);
  return getSpline(knot, d[0], d[1], d[2], d[3], u);
}

template <typename T>
GPUdi() GPUgeneric() T CompactSplineIrregular1D::getSpline(GPUgeneric() const T data[], float u) const
{
  /// Get interpolated value for f(u) using data array data[2*getNumberOfKnots()] = {F0,D0,..,Fn-1,Dn-1}
  int iknot = getKnotIndex(u);
  const CompactSplineIrregular1D::Knot& knot = getKnot(iknot);
  const T* d = data + (iknot + iknot);
  return getSpline(knot, d[0], d[1], d[2], d[3], u);
}

GPUdi() int CompactSplineIrregular1D::getKnotIndex(float u) const
{
  /// Get i: u is in [knot_i, knot_{i+1}) interval
  /// no border check! u must be in [0,mUmax]
  return getBin2KnotMap()[(int)u];
}

GPUdi() int CompactSplineIrregular1D::getKnotIndexSafe(float u) const
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

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif