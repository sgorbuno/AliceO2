// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  CompactSplineIrregular2D3D.cxx
/// \brief Implementation of CompactSplineIrregular2D3D class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#include "CompactSplineIrregular2D3D.h"

#if !defined(GPUCA_GPUCODE)
#include <iostream>
#endif

using namespace GPUCA_NAMESPACE::gpu;

CompactSplineIrregular2D3D::CompactSplineIrregular2D3D() : FlatObject(), mGridU(), mGridV()
{
  /// Default constructor. Creates an empty uninitialised object
}

void CompactSplineIrregular2D3D::destroy()
{
  /// See FlatObject for description
  mGridU.destroy();
  mGridV.destroy();
  FlatObject::destroy();
}

void CompactSplineIrregular2D3D::cloneFromObject(const CompactSplineIrregular2D3D& obj, char* newFlatBufferPtr)
{
  /// See FlatObject for description

  const char* oldFlatBufferPtr = obj.mFlatBufferPtr;

  FlatObject::cloneFromObject(obj, newFlatBufferPtr);

  char* bufferU = FlatObject::relocatePointer(oldFlatBufferPtr, mFlatBufferPtr, obj.mGridU.getFlatBufferPtr());
  char* bufferV = FlatObject::relocatePointer(oldFlatBufferPtr, mFlatBufferPtr, obj.mGridV.getFlatBufferPtr());

  mGridU.cloneFromObject(obj.mGridU, bufferU);
  mGridV.cloneFromObject(obj.mGridV, bufferV);
}

void CompactSplineIrregular2D3D::moveBufferTo(char* newFlatBufferPtr)
{
  /// See FlatObject for description
#ifndef GPUCA_GPUCODE
  char* oldFlatBufferPtr = mFlatBufferPtr;
  FlatObject::moveBufferTo(newFlatBufferPtr);
  char* currFlatBufferPtr = mFlatBufferPtr;
  mFlatBufferPtr = oldFlatBufferPtr;
  setActualBufferAddress(currFlatBufferPtr);
#endif
}

void CompactSplineIrregular2D3D::setActualBufferAddress(char* actualFlatBufferPtr)
{
  /// See FlatObject for description
  char* bufferU = FlatObject::relocatePointer(mFlatBufferPtr, actualFlatBufferPtr, mGridU.getFlatBufferPtr());
  char* bufferV = FlatObject::relocatePointer(mFlatBufferPtr, actualFlatBufferPtr, mGridV.getFlatBufferPtr());
  mGridU.setActualBufferAddress(bufferU);
  mGridV.setActualBufferAddress(bufferV);
  FlatObject::setActualBufferAddress(actualFlatBufferPtr);
}

void CompactSplineIrregular2D3D::setFutureBufferAddress(char* futureFlatBufferPtr)
{
  /// See FlatObject for description
  char* bufferU = relocatePointer(mFlatBufferPtr, futureFlatBufferPtr, mGridU.getFlatBufferPtr());
  char* bufferV = relocatePointer(mFlatBufferPtr, futureFlatBufferPtr, mGridV.getFlatBufferPtr());
  mGridU.setFutureBufferAddress(bufferU);
  mGridV.setFutureBufferAddress(bufferV);
  FlatObject::setFutureBufferAddress(futureFlatBufferPtr);
}

void CompactSplineIrregular2D3D::construct(int numberOfKnotsU, const float knotsU[], int numberOfAxisBinsU, int numberOfKnotsV, const float knotsV[], int numberOfAxisBinsV)
{
  /// Constructor
  ///
  /// Number of knots created and their values may differ from the input values:
  /// - Edge knots 0.f and 1.f will be added if they are not present.
  /// - Knot values are rounded to closest axis bins: k*1./numberOfAxisBins.
  /// - Knots rounded to the same axis bin will be merged
  /// - At least 2 knots and at least 1 axis bin will be created in both directions
  ///
  /// \param numberOfKnotsU     U axis: Number of knots in knots[] array
  /// \param knotsU             U axis: Array of knots.
  /// \param numberOfAxisBinsU  U axis: Number of axis bins to map arbitrary U coordinate to
  ///                           an appropriate [knot(i),knot(i+1)] interval.
  ///                           The knot positions have "granularity" of 1./numberOfAxisBins
  ///
  /// \param numberOfKnotsV     V axis: Number of knots in knots[] array
  /// \param knotsV             V axis: Array of knots.
  /// \param numberOfAxisBinsV  V axis: Number of axis bins to map U coordinate to
  ///                           an appropriate [knot(i),knot(i+1)] interval.
  ///                           The knot positions have "granularity" of 1./numberOfAxisBins
  ///

  FlatObject::startConstruction();

  //mGridU.construct(numberOfKnotsU, knotsU, numberOfAxisBinsU);
  //mGridV.construct(numberOfKnotsV, knotsV, numberOfAxisBinsV);

  size_t vOffset = alignSize(mGridU.getFlatBufferSize(), mGridV.getBufferAlignmentBytes());

  FlatObject::finishConstruction(vOffset + mGridV.getFlatBufferSize());

  mGridU.moveBufferTo(mFlatBufferPtr);
  mGridV.moveBufferTo(mFlatBufferPtr + vOffset);
}

void CompactSplineIrregular2D3D::constructRegular(int numberOfKnotsU, int numberOfKnotsV)
{
  /// Constructor for a regular spline
  /// \param numberOfKnotsU     U axis: Number of knots in knots[] array
  /// \param numberOfKnotsV     V axis: Number of knots in knots[] array
  ///

  FlatObject::startConstruction();

  mGridU.constructRegular(numberOfKnotsU);
  mGridV.constructRegular(numberOfKnotsV);

  size_t vOffset = alignSize(mGridU.getFlatBufferSize(), mGridV.getBufferAlignmentBytes());

  FlatObject::finishConstruction(vOffset + mGridV.getFlatBufferSize());

  mGridU.moveBufferTo(mFlatBufferPtr);
  mGridV.moveBufferTo(mFlatBufferPtr + vOffset);
}

void CompactSplineIrregular2D3D::print() const
{
#if !defined(GPUCA_GPUCODE)
  std::cout << " Irregular Spline 2D3D: " << std::endl;
  std::cout << " grid U: " << std::endl;
  mGridU.print();
  std::cout << " grid V: " << std::endl;
  mGridV.print();
#endif
}
