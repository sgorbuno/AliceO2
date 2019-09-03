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

#ifndef __OPENCL__
#include <cstddef>
#include <memory>
#include <cstring>
#endif

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
/// The interpolation is performed on n knots {U0==0., U1, .., Un<1.}
/// with given function values {F0, ..., Fn} and derivatives [Z0,..,Zn] at the knots.
///
/// An interpolation in each interval between two knots is performed by 3-th degree polynom.
/// The polynoms cross the Fi values at the knots and have contnious 1-st derivative.
/// Depending on the initialization of the first derivatives, the second derivative may not be continious.
///
/// The knots should belong to half-interval [0,1), the distance between knots is (almost) arbitrary.
///
/// Nothing which depends on F is stored in the class.
/// Therefore one can use the same class for different F functions.
/// The function values {F0,..,Fn} and the derivatives {Z0,..,Zn} have to be provided by user for each call.
///
/// The class performs a fast search of a spline interval: (float U ) -> [int iKnot, int iKnot+1 ).
/// For that purpose, initial U coordinates of the knots are rounded to the closest i*1./nAxisBins values.
/// Number of knots and they U coordinates may change during initialisation!
///
/// The minimal number of knots is 2, the minimal number of axis bins is 1
///
/// Knots U0=0. and Un=1. are always present. They are added automatically when they are not set by an user.
///
/// User should provide function values Fi and the derivatives Zi for all !constructed! knots.
///
/// ------------
///
///
///  Example of creating a spline:
///
///  const int nKnots=2;
///  float knots[nKnots] = {0., 1.};
///  CompactSplineIrregular1D spline;
///  spline.construct(nKnots, knots, 1);
///  float data[2*nKnots] = { 3.5, 0.01, 2.0, -0.01};
///  spline.getSpline( f, 0.0 ); // == 3.5
///  spline.getSpline( f, 0.1 ); // == some interpolated value
///  spline.getSpline( f, 0.5 ); // == some interpolated value
///  spline.getSpline( f, 1.0 ); // == 2.0
///
class CompactSplineIrregular1D : public FlatObject
{
 public:
  ///
  /// \brief The struct represents a knot(i) and interval [ knot(i), knot(i+1) ]
  ///
  struct Knot {
    float u;  ///< u coordinate of the knot i
    float L;  ///< length of the [i, i+1] segment
    float Li; ///< inverse length of the [i, i+1] segment
  };

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
#ifndef GPUCA_GPUCODE
  using FlatObject::moveBufferTo;
#endif

  /// Moving the class with its external buffer to another location

  using FlatObject::setActualBufferAddress;
  using FlatObject::setFutureBufferAddress;

  /// _______________  Construction interface  ________________________

  /// Constructor
  ///
  /// Number of knots created and their values may differ from the input values:
  /// - Edge knots 0.f and 1.f will be added if they are not present.
  /// - Knot values are rounded to closest axis bins: k*1./numberOfAxisBins.
  /// - Knots which are too close to each other will be merged
  /// - At least 5 knots and at least 4 axis bins will be created for consistency reason
  ///
  /// \param numberOfKnots     Number of knots in knots[] array
  /// \param knots             Array of knots.
  /// \param numberOfAxisBins Number of axis bins to map U coordinate to
  ///                          an appropriate [knot(i),knot(i+1)] interval.
  ///                          The knot positions have a "granularity" of 1./numberOfAxisBins
  ///
  void construct(int numberOfKnots, const float knots[], int numberOfAxisBins);

  /// Constructor for a regular spline
  void constructRegular(int numberOfKnots);

  /// _______________  Main functionality   ________________________

  /// Get interpolated value for f(u) using spline at the interval [knot,next knot] with function values f0, f1 and derivatives z0, z1
  template <typename T>
  GPUd() static T getSpline(const CompactSplineIrregular1D::Knot& knot, T f0, T z0, T f1, T z1, float u);

  /// Get interpolated value for f(u) using data array data[2*getNumberOfKnots()]
  template <typename T>
  GPUd() T getSpline(const T data[], float u) const;

  /// Get number of knots
  GPUd() int getNumberOfKnots() const { return mNumberOfKnots; }

  /// Get index of associated knot for a given U coordinate.
  ///
  /// Note: U values from the first interval are mapped to the second inrerval.
  /// Values from the last interval are mapped to the previous interval.
  ///
  GPUd() int getKnotIndex(float u) const;

  /// Get i-th knot, no border check performed!
  GPUd() const CompactSplineIrregular1D::Knot& getKnot(int i) const { return getKnots()[i]; }

  /// Get array of knots
  GPUd() const CompactSplineIrregular1D::Knot* getKnots() const { return reinterpret_cast<const CompactSplineIrregular1D::Knot*>(mFlatBufferPtr); }

  /// Get minimal required alignment for the class
  static constexpr size_t getClassAlignmentBytes() { return 8; }

  /// Get minimal required alignment for the flat buffer
  static constexpr size_t getBufferAlignmentBytes() { return 8; }

  /// Get minimal required alignment for the spline data
  static constexpr size_t getDataAlignmentBytes() { return 8; }

  /// technical stuff

  /// Get a map  (U axis bin index) -> (corresponding knot index)
  GPUd() const int* getBin2KnotMap() const { return reinterpret_cast<const int*>(mFlatBufferPtr + mBin2KnotMapOffset); }

  /// Get number of axis bins
  int getNumberOfAxisBins() const { return mNumberOfAxisBins; }

  /// Print method
  void print() const;

 private:
  /// Non-const accessor to knots array
  CompactSplineIrregular1D::Knot* getKnotsNonConst() { return reinterpret_cast<CompactSplineIrregular1D::Knot*>(mFlatBufferPtr); }

  /// Non-const accessor to bins->knots map
  int* getBin2KnotMapNonConst() { return reinterpret_cast<int*>(mFlatBufferPtr + mBin2KnotMapOffset); }

  ///
  /// ====  Data members   ====
  ///

  int mNumberOfKnots;              ///< n knots on the grid
  int mNumberOfAxisBins;           ///< number of axis bins
  unsigned int mBin2KnotMapOffset; ///< pointer to (axis bin) -> (knot) map in mFlatBufferPtr array
};

/// ====================================================
///       Inline implementations of some methods
/// ====================================================

template <typename T>
GPUdi() T CompactSplineIrregular1D::getSpline(const CompactSplineIrregular1D::Knot& knot0, T f0, T z0, T f1, T z1, float u)
{
  /// static method
  /// Get interpolated value for f(u) using spline at knot "knot1" and function values at knots {knot_0,knot_1,knot_2,knot_3}

  T x = T((u - knot0.u) * knot0.Li); // scaled u
  T x2 = x * x;
  /*
  T xm1 = x-1;
  float cf1 = x2*(3-2*x);
  float cf0 = 1-cf1;
  float cz0 = x*xm1*xm1*knot0.L;
  float cz1 = x2*xm1*knot0.L;
  return cf0*f0 + cf1*f1 + cz0*z0 + cz1*z1;
  */

  z0 *= knot0.L; // scaled u derivative at the knot 0
  z1 *= knot0.L; // scaled u derivative at the knot 1
  f1 -= f0;

  // f(x) = ax^3 + bx^2 + cx + d
  //
  // f(0) = f0
  // f(1) = f1
  // f'(0) = z0
  // f'(1) = z1
  //

  // T d = f0;
  // T c = z0;
  T a = -f1 - f1 + z0 + z1;
  T b = f1 - z0 - a;
  return a * x * x2 + b * x2 + z0 * x + f0;
}

template <typename T>
GPUdi() T CompactSplineIrregular1D::getSpline(const T data[], float u) const
{
  /// Get interpolated value for f(u) using data array data[2*getNumberOfKnots()]
  int iknot = getKnotIndex(u);
  const CompactSplineIrregular1D::Knot& knot = getKnot(iknot);
  const T* d = data + 2 * iknot;
  return getSpline(knot, d[0], d[1], d[2], d[3], u);
}

GPUdi() int CompactSplineIrregular1D::getKnotIndex(float u) const
{
  /// get i: u is in [knot_i, knot_{i+1})
  int ibin = (int)(u * mNumberOfAxisBins);
  if (ibin < 0) {
    ibin = 0;
  }
  if (ibin > mNumberOfAxisBins - 1) {
    ibin = mNumberOfAxisBins - 1;
  }
  return getBin2KnotMap()[ibin];
}

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
