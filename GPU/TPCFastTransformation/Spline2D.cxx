// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  CompactSpline2D.cxx
/// \brief Implementation of CompactSpline2D class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#include "CompactSpline2D.h"

#if !defined(GPUCA_GPUCODE)
#include <iostream>
#endif

using namespace GPUCA_NAMESPACE::gpu;

CompactSpline2D::CompactSpline2D() : FlatObject(), mGridU(), mGridV()
{
  /// Default constructor. Creates an empty uninitialised object
}

void CompactSpline2D::destroy()
{
  /// See FlatObject for description
  mGridU.destroy();
  mGridV.destroy();
  FlatObject::destroy();
}

void CompactSpline2D::cloneFromObject(const CompactSpline2D& obj, char* newFlatBufferPtr)
{
  /// See FlatObject for description

  const char* oldFlatBufferPtr = obj.mFlatBufferPtr;

  FlatObject::cloneFromObject(obj, newFlatBufferPtr);

  char* bufferU = FlatObject::relocatePointer(oldFlatBufferPtr, mFlatBufferPtr, obj.mGridU.getFlatBufferPtr());
  char* bufferV = FlatObject::relocatePointer(oldFlatBufferPtr, mFlatBufferPtr, obj.mGridV.getFlatBufferPtr());

  mGridU.cloneFromObject(obj.mGridU, bufferU);
  mGridV.cloneFromObject(obj.mGridV, bufferV);
}

void CompactSpline2D::moveBufferTo(char* newFlatBufferPtr)
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

void CompactSpline2D::setActualBufferAddress(char* actualFlatBufferPtr)
{
  /// See FlatObject for description
  char* bufferU = FlatObject::relocatePointer(mFlatBufferPtr, actualFlatBufferPtr, mGridU.getFlatBufferPtr());
  char* bufferV = FlatObject::relocatePointer(mFlatBufferPtr, actualFlatBufferPtr, mGridV.getFlatBufferPtr());
  mGridU.setActualBufferAddress(bufferU);
  mGridV.setActualBufferAddress(bufferV);
  FlatObject::setActualBufferAddress(actualFlatBufferPtr);
}

void CompactSpline2D::setFutureBufferAddress(char* futureFlatBufferPtr)
{
  /// See FlatObject for description
  char* bufferU = relocatePointer(mFlatBufferPtr, futureFlatBufferPtr, mGridU.getFlatBufferPtr());
  char* bufferV = relocatePointer(mFlatBufferPtr, futureFlatBufferPtr, mGridV.getFlatBufferPtr());
  mGridU.setFutureBufferAddress(bufferU);
  mGridV.setFutureBufferAddress(bufferV);
  FlatObject::setFutureBufferAddress(futureFlatBufferPtr);
}

void CompactSpline2D::construct(int numberOfKnotsU, const int knotsU[], int numberOfKnotsV, const int knotsV[])
{
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

  FlatObject::startConstruction();

  mGridU.construct(numberOfKnotsU, knotsU);
  mGridV.construct(numberOfKnotsV, knotsV);

  size_t vOffset = alignSize(mGridU.getFlatBufferSize(), mGridV.getBufferAlignmentBytes());

  FlatObject::finishConstruction(vOffset + mGridV.getFlatBufferSize());

  mGridU.moveBufferTo(mFlatBufferPtr);
  mGridV.moveBufferTo(mFlatBufferPtr + vOffset);
}

void CompactSpline2D::constructRegular(int numberOfKnotsU, int numberOfKnotsV)
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

void CompactSpline2D::print() const
{
#if !defined(GPUCA_GPUCODE)
  std::cout << " Irregular Spline 2D: " << std::endl;
  std::cout << " grid U: " << std::endl;
  mGridU.print();
  std::cout << " grid V: " << std::endl;
  mGridV.print();
#endif
}
