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
/// The spline interpolates a generic function F:[0,1]->R.
///
/// Let's call the function parameter U, the function value F.
/// The interpolation is performed on n knots {U0==0., U1, .., Un-1==1.}
/// with given function values {F0, ..., Fn-1} and derivatives [Z0,..,Zn-1] at the knots.
///
/// An interpolation on each segment between two knots is performed by a 3-th degree polynom.
/// The polynoms and they 1-st derivatives are continuous at the knots.
/// Depending on the initialization of the first derivative, the second derivative may or may not be continious.
///
/// The knots should belong to a segment [0,1], the distance between knots is (almost) arbitrary.
///
/// Nothing which depends on F is stored in the class,
/// therefore one can use the same class for interpolation of different input functions.
/// The function values {F0,..,Fn-1} and the derivatives {Z0,..,Zn-1} have to be provided by user for each call.
///
/// An utility for searching of the spline segment: (float U ) -> [int iKnot, int iKnot+1 ] is implemented in a fast way.
/// For that purpose, initial U coordinates of the knots are rounded to the nearby i*1./nAxisBins values.
/// Due to this implementation, the number of knots and they U coordinates may change during initialization!
///
/// The minimal number of knots is 2, the minimal number of axis bins is 1
///
/// Knots U0=0. and Un=1. are always present. They are added automatically when they are not set by the user.
///
/// The user should provide function values Fi and the derivatives Zi for all !constructed! knots.
/// They can be calculated using corresponding utilities in CompactSplineHelper class.
///
/// ------------
///
///  Example of creating a spline:
///
///  const int nKnots=3;
///  float knots[nKnots] = {0., 0.2., 1.}; // original knot positions
///  const int nAxisBins = 3;
///  CompactSplineIrregular1D spline;
///  spline.construct(nKnots, knots, nAxisBins); // knot positions will be rounded to 1/nAxisBins:  {0, 0.3, 1.}
///  float data[2*nKnots] = { 3.5, 0.01, 2.0, -0.01, 3.1, 0.02};
///  spline.getSpline( f, 0.0 ); // == 3.5
///  spline.getSpline( f, 0.1 ); // == some interpolated value
///  spline.getSpline( f, 0.3 ); // == 2.0
///  spline.getSpline( f, 1.0 ); // == 3.1
///
class CompactSplineIrregular1D : public FlatObject
{
 public:
  ///
  /// \brief The struct represents a knot(i) and interval [ knot(i), knot(i+1) ]
  ///
  struct Knot {
    float u;  ///< u coordinate of the knot i
    float Li; ///< inverse length of the [i, i+1] segment
  };

  /// _____________  Version control __________________________

  static constexpr int getVersion() { return 1; }
  size_t getDataSize() const { return 2 * mNumberOfKnots * sizeof(float); }

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

  void cloneFromObject(const CompactSplineIrregular1D& obj, char* newFlatBufferPtr);
  void destroy();

  /// Making the data buffer external

  using FlatObject::releaseInternalBuffer;
  void moveBufferTo(char* newBufferPtr);

  /// Moving the class with its external buffer to another location

  void setActualBufferAddress(char* actualFlatBufferPtr);
  void setFutureBufferAddress(char* futureFlatBufferPtr);

  /// Get minimal required alignment for the spline data
  static constexpr size_t getDataAlignmentBytes() { return 2 * sizeof(float); }

  /// _______________  Construction interface  ________________________

  /// Constructor
  ///
  /// Number of created knots and their positions may differ from the input values:
  /// - Edge knots 0.f and 1.f will be added if they are not present.
  /// - Knot positions are rounded to the nearby axis bins: k*1./numberOfAxisBins.
  /// - Knots rounded to the same axis bin will be merged
  /// - At least 2 knots and at least 1 axis bin will be created
  ///
  /// \param numberOfKnots     Number of knots in knots[] array
  /// \param knots             Array of knots.
  /// \param numberOfAxisBins Number of axis bins to map U coordinate to
  ///                          an appropriate [knot(i),knot(i+1)] interval.
  ///                          The knot positions have a "granularity" of 1./numberOfAxisBins
  ///
  void construct(int numberOfKnots, const float knots[], int numberOfAxisBins);

  /// Constructor for a regular spline. Knots will be placed at positions i/(numberOfKnots-1)
  void constructRegular(int numberOfKnots);

  /// _______________  Main functionality   ________________________

  /// Get interpolated value for F(u) using spline at the segment [knot,next knot] with function values f0, f1 and derivatives z0, z1
  template <typename T>
  GPUd() static T getSpline(const CompactSplineIrregular1D::Knot& knot, T f0, T z0, T f1, T z1, float u);

  /// Get interpolated value for F(u) using data array data[2*getNumberOfKnots()] == {F0,Z0, .., Fn-1,Zn-1}
  template <typename T>
  GPUd() T getSpline(const T data[], float u) const;

  /// Get interpolated value for F(u) using data array data[2*getNumberOfKnots()] == {F0,Z0, .., Fn-1,Zn-1}, with a border check
  template <typename T>
  GPUd() T getSplineSafe(const T data[], float u) const;

  /// Get number of knots
  GPUd() int getNumberOfKnots() const { return mNumberOfKnots; }

  /// Get index of associated knot for a given U coordinate.
  GPUd() int getKnotIndex(float u) const;

  /// Get i-th knot, no border check performed!
  GPUd() const CompactSplineIrregular1D::Knot& getKnot(int i) const { return getKnots()[i]; }

  /// Get the array of knots
  GPUd() const CompactSplineIrregular1D::Knot* getKnots() const { return reinterpret_cast<const CompactSplineIrregular1D::Knot*>(mFlatBufferPtr); }

  /// technical stuff

  /// Get a map  (U axis bin index) -> (corresponding knot index)
  GPUd() const int* getBin2KnotMap() const { return mBin2KnotMap; }

  /// Get number of axis bins
  int getNumberOfAxisBins() const { return mNumberOfAxisBins; }

  /// Print method
  void print() const;

 private:
  /// Non-const accessor to knots array
  CompactSplineIrregular1D::Knot* getKnotsNonConst() { return reinterpret_cast<CompactSplineIrregular1D::Knot*>(mFlatBufferPtr); }

  /// Non-const accessor to bins->knots map
  int* getBin2KnotMapNonConst() { return mBin2KnotMap; }

  ///
  /// ====  Data members   ====
  ///

  int mNumberOfKnots;    ///< n knots on the grid
  int mNumberOfAxisBins; ///< number of axis bins
  int* mBin2KnotMap;     ///< pointer to (axis bin) -> (knot) map inside the mFlatBufferPtr array
};

/// ====================================================
///       Inline implementations of some methods
/// ====================================================

template <typename T>
GPUdi() T CompactSplineIrregular1D::getSpline(const CompactSplineIrregular1D::Knot& knot0, T f0, T z0, T f1, T z1, float u)
{
  /// static method
  /// Get interpolated value for f(u) using spline at knot "knot0" and function & derivative values at knots {knot_0,knot_1}

  T uu = T(u - knot0.u);
  T x = uu * T(knot0.Li); // scaled u

  /* another way to calculate
  T xm1 = x-1;
  T x2 = x * x;
  float cf1 = x2*(3-2*x);
  float cf0 = 1-cf1;
  float cz0 = x*xm1*xm1*knot0.L;
  float cz1 = x2*xm1*knot0.L;
  return cf0*f0 + cf1*f1 + cz0*z0 + cz1*z1;
  */
  f1 = (f1 - f0) * knot0.Li;
  T a = -f1 - f1 + z0 + z1;
  T b = f1 - z0 - a;
  return ((a * x + b) * x + z0) * uu + f0;
}

template <typename T>
GPUdi() T CompactSplineIrregular1D::getSplineSafe(const T data[], float u) const
{
  /// Get interpolated value for f(u) using data array data[2*getNumberOfKnots()]
  if (u < 0.f)
    u = 0.f;
  else if (u > 1.f)
    u = 1.f;
  return getSpline(data, u);
}

template <typename T>
GPUdi() T CompactSplineIrregular1D::getSpline(const T data[], float u) const
{
  /// get interpolated value for f(u) using data array data[2*getNumberOfKnots()] = {F0,Z0,..,Fn-1,Zn-1}
  int iknot = getKnotIndex(u);
  const CompactSplineIrregular1D::Knot& knot = getKnot(iknot);
  const T* d = data + (iknot + iknot);
  return getSpline(knot, d[0], d[1], d[2], d[3], u);
}

GPUdi() int CompactSplineIrregular1D::getKnotIndex(float u) const
{
  /// get i: u is in [knot_i, knot_{i+1})
  int ibin = (int)(u * mNumberOfAxisBins);
  return getBin2KnotMap()[ibin];
}

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif