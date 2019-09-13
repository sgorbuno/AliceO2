// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  CompactSplineIrregular2D3D.h
/// \brief Definition of CompactSplineIrregular2D3D class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#ifndef ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_COMPACTSPLINEIRREGULAR2D3D_H
#define ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_COMPACTSPLINEIRREGULAR2D3D_H

#include "CompactSplineIrregular1D.h"
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
/// The CompactSplineIrregular2D3D class represents twoo-dimensional spline interpolation on nonunifom (irregular) grid.
///
/// The class is flat C structure. No virtual methods, no ROOT types are used.
/// It is designed for spline parameterisation of TPC transformation.
///
/// ---
/// The spline interpolates a generic function F:[u,v)->(x,y,z),
/// where u,v belong to [0,1]x[0,1]
///
/// It is an extension of CompactSplineIrregular1D class, see CompactSplineIrregular1D.h for more details.
///
/// Important:
///   -- The number of knots and their positions may change during initialisation
///
/// ------------
///
///  Example of creating a spline:
///
///  const int nKnotsU=5;
///  const int nKnotsV=6;
///  float knotsU[nKnotsU] = {0., 0.25, 0.5, 0.7, 1.};
///  float knotsV[nKnotsV] = {0., 0.2, 0.3, 0.45, 0.8, 1.};
///  CompactSplineIrregular2D3D spline(nKnotsU, knotsU, 4, nKnotsV, knotsV, 10 );
///  float f[nKnotsU*nKnotsV*3] = { 3.5, .... };
///  spline.getSpline( f, 0., 0. ); // == 3.5
///  spline.getSpline( f, 0.1, 0.32 ); // == some interpolated value
///
class CompactSplineIrregular2D3D : public FlatObject
{
 public:
  /// _____________  Version control __________________________

  /// Version number
  GPUd() static constexpr int getVersion() { return 1; }

  /// Size of the data array in elements, must be multiplied by sizeof(float)
  GPUd() size_t getDataSizeInelements() const { return 4 * getNumberOfKnots(); }

  /// _____________  Constructors / destructors __________________________

  /// Default constructor. Creates an empty uninitialised object
  CompactSplineIrregular2D3D();

  /// Copy constructor: disabled to avoid ambiguity. Use cloneFromObject() instead
  CompactSplineIrregular2D3D(const CompactSplineIrregular2D3D&) CON_DELETE;

  /// Assignment operator: disabled to avoid ambiguity. Use cloneFromObject() instead
  CompactSplineIrregular2D3D& operator=(const CompactSplineIrregular2D3D&) CON_DELETE;

  /// Destructor
  ~CompactSplineIrregular2D3D() CON_DEFAULT;

  /// _____________  FlatObject functionality, see FlatObject class for description  ____________

  /// Memory alignment

  using FlatObject::getBufferAlignmentBytes;
  using FlatObject::getClassAlignmentBytes;
  /// Get minimal required alignment for the spline data
  static constexpr size_t getDataAlignmentBytes() { return 4 * sizeof(float); }

  /// Construction interface

  void cloneFromObject(const CompactSplineIrregular2D3D& obj, char* newFlatBufferPtr);
  void destroy();

  /// Making the data buffer external

  using FlatObject::releaseInternalBuffer;
  void moveBufferTo(char* newBufferPtr);

  /// Moving the class with its external buffer to another location

  void setActualBufferAddress(char* actualFlatBufferPtr);
  void setFutureBufferAddress(char* futureFlatBufferPtr);

  /// _______________  Construction interface  ________________________

  /// Constructor
  ///
  /// Number of knots created and their values may differ from the input values:
  /// - Edge knots 0.f and 1.f will be added if they are not present.
  /// - Knot values are rounded to closest axis bins: k*1./numberOfAxisBins.
  /// - Knots which are too close to each other will be merged
  /// - At least 2 knots and at least 1 axis bin will be created in both directions
  ///
  /// \param numberOfKnotsU     U axis: Number of knots in knots[] array
  /// \param knotsU             U axis: Array of knots.
  /// \param numberOfAxisBinsU  U axis: Number of axis bins to map U coordinate to
  ///                           an appropriate [knot(i),knot(i+1)] interval.
  ///                           The knot positions have a "granularity" of 1./numberOfAxisBins
  ///
  /// \param numberOfKnotsV     V axis: Number of knots in knots[] array
  /// \param knotsV             V axis: Array of knots.
  /// \param numberOfAxisBinsV  V axis: Number of axis bins to map U coordinate to
  ///                           an appropriate [knot(i),knot(i+1)] interval.
  ///                           The knot positions have a "granularity" of 1./numberOfAxisBins
  ///
  void construct(int numberOfKnotsU, const float knotsU[], int numberOfAxisBinsU, int numberOfKnotsV, const float knotsV[], int numberOfAxisBinsV);

  /// Constructor for a regular spline
  void constructRegular(int numberOfKnotsU, int numberOfKnotsV);

  /// _______________  Main functionality   ________________________

  /// Get interpolated value for f(u,v) using data array data[4*getNumberOfKnots()] with corrected edges
  template <typename T>
  GPUd() void getSpline(GPUgeneric() const T* data, float u, float v, GPUgeneric() T& x, GPUgeneric() T& y, GPUgeneric() T& z) const;

  /// Same as getSpline, but using vectorized calculation.
  /// \param correctedData should be at least 128-bit aligned
  GPUd() void getSplineVec(const float* data, float u, float v, float& x, float& y, float& z) const;

  /// Get number total of knots: UxV
  GPUd() int getNumberOfKnots() const { return mGridU.getNumberOfKnots() * mGridV.getNumberOfKnots(); }

  /// Get 1-D grid for U coordinate
  GPUd() const CompactSplineIrregular1D& getGridU() const { return mGridU; }

  /// Get 1-D grid for V coordinate
  GPUd() const CompactSplineIrregular1D& getGridV() const { return mGridV; }

  /// Get 1-D grid for U or V coordinate
  GPUd() const CompactSplineIrregular1D& getGrid(int uv) const { return (uv == 0) ? mGridU : mGridV; }

  /// Get u,v of i-th knot
  GPUd() void getKnotUV(int iKnot, float& u, float& v) const;

  /// technical stuff

  /// Get offset of GridU flat data in the flat buffer
  size_t getGridUOffset() const { return mGridU.getFlatBufferPtr() - mFlatBufferPtr; }

  /// Get offset of GridV flat data in the flat buffer
  size_t getGridVOffset() const { return mGridV.getFlatBufferPtr() - mFlatBufferPtr; }

  /// Print method
  void print() const;

 private:
  ///
  /// ====  Data members   ====
  ///

  CompactSplineIrregular1D mGridU; ///< grid for U axis
  CompactSplineIrregular1D mGridV; ///< grid for V axis
};

/// ====================================================
///       Inline implementations of some methods
/// ====================================================

GPUdi() void CompactSplineIrregular2D3D::getKnotUV(int iKnot, float& u, float& v) const
{
  /// Get u,v of i-th knot
  const CompactSplineIrregular1D& gridU = getGridU();
  const CompactSplineIrregular1D& gridV = getGridV();
  int nu = gridU.getNumberOfKnots();
  int iv = iKnot / nu;
  int iu = iKnot % nu;
  u = gridU.getKnot(iu).u;
  v = gridV.getKnot(iv).u;
}

template <typename T>
GPUdi() void CompactSplineIrregular2D3D::getSpline(GPUgeneric() const T* correctedData, float u, float v, GPUgeneric() T& x, GPUgeneric() T& y, GPUgeneric() T& z) const
{
  // Get interpolated value for f(u,v) using data array correctedData[getNumberOfKnots()] with corrected edges

  const CompactSplineIrregular1D& gridU = getGridU();
  const CompactSplineIrregular1D& gridV = getGridV();
  int nu = gridU.getNumberOfKnots();
  int iu = gridU.getKnotIndex(u);
  int iv = gridV.getKnotIndex(v);

  const CompactSplineIrregular1D::Knot& knotU = gridU.getKnot(iu);
  const CompactSplineIrregular1D::Knot& knotV = gridV.getKnot(iv);

  const T* dataV0 = correctedData + (nu * (iv - 1) + iu - 1) * 3;
  const T* dataV1 = dataV0 + 3 * nu;
  const T* dataV2 = dataV0 + 6 * nu;
  const T* dataV3 = dataV0 + 9 * nu;

  T dataV[12];
  for (int i = 0; i < 12; i++) {
    dataV[i] = gridV.getSpline(knotV, dataV0[i], dataV1[i], dataV2[i], dataV3[i], v);
  }

  T* dataU0 = dataV + 0;
  T* dataU1 = dataV + 3;
  T* dataU2 = dataV + 6;
  T* dataU3 = dataV + 9;

  T res[3];
  for (int i = 0; i < 3; i++) {
    res[i] = gridU.getSpline(knotU, dataU0[i], dataU1[i], dataU2[i], dataU3[i], u);
  }
  x = res[0];
  y = res[1];
  z = res[2];
}

GPUdi() void CompactSplineIrregular2D3D::getSplineVec(const float* correctedData, float u, float v, float& x, float& y, float& z) const
{
// Same as getSpline, but using vectorized calculation.
// \param correctedData should be at least 128-bit aligned

#if !defined(__CINT__) && !defined(__ROOTCINT__) && !defined(GPUCA_GPUCODE) && !defined(GPUCA_NO_VC) && defined(__cplusplus) && __cplusplus >= 201703L
  const CompactSplineIrregular1D& gridU = getGridU();
  const CompactSplineIrregular1D& gridV = getGridV();
  int nu = gridU.getNumberOfKnots();
  int iu = gridU.getKnotIndex(u);
  int iv = gridV.getKnotIndex(v);
  const CompactSplineIrregular1D::Knot& knotU = gridU.getKnot(iu);
  const CompactSplineIrregular1D::Knot& knotV = gridV.getKnot(iv);

  const float* dataV0 = correctedData + (nu * (iv - 1) + iu - 1) * 3;
  const float* dataV1 = dataV0 + 3 * nu;
  const float* dataV2 = dataV0 + 6 * nu;
  const float* dataV3 = dataV0 + 9 * nu;

  // calculate F values at V==v and U == Ui of knots

  Vc::SimdArray<float, 12> dataV0vec(dataV0), dataV1vec(dataV1), dataV2vec(dataV2), dataV3vec(dataV3), dataVvec = gridV.getSpline(knotV, dataV0vec, dataV1vec, dataV2vec, dataV3vec, v);

  using V = std::conditional_t<(Vc::float_v::size() >= 4),
                               Vc::float_v,
                               Vc::SimdArray<float, 3>>;

  float dataV[9 + V::size()];
  // dataVvec.scatter( dataV, Vc::SimdArray<float,12>::IndexType::IndexesFromZero() );
  //dataVvec.scatter(dataV, Vc::SimdArray<uint, 12>(Vc::IndexesFromZero));
  dataVvec.store(dataV, Vc::Unaligned);

  for (unsigned int i = 12; i < 9 + V::size(); i++) // fill not used part of the vector with 0
    dataV[i] = 0.f;

  // calculate F values at V==v and U == u
  V dataU0vec(dataV), dataU1vec(dataV + 3), dataU2vec(dataV + 6), dataU3vec(dataV + 9);

  V dataUVvec = gridU.getSpline(knotU, dataU0vec, dataU1vec, dataU2vec, dataU3vec, u);

  x = dataUVvec[0];
  y = dataUVvec[1];
  z = dataUVvec[2];
#else
  getSpline(correctedData, u, v, x, y, z);
#endif
}
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
