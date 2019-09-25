// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  CompactSplineIrregular1D.cxx
/// \brief Implementation of CompactSplineIrregular1D class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#include "CompactSplineIrregular1D.h"
#include <cmath>

#if !defined(GPUCA_GPUCODE)
#include <vector>
#include <algorithm>
#include <iostream>
#endif

using namespace GPUCA_NAMESPACE::gpu;

CompactSplineIrregular1D::CompactSplineIrregular1D() : FlatObject(), mNumberOfKnots(0), mUmax(0), mBin2KnotMap(0)
{
  /// Default constructor. Creates an empty uninitialised object
}

void CompactSplineIrregular1D::destroy()
{
  /// See FlatObject for description
  mNumberOfKnots = 0;
  mUmax = 0;
  mBin2KnotMap = nullptr;
  FlatObject::destroy();
}

#if !defined(GPUCA_GPUCODE)
void CompactSplineIrregular1D::cloneFromObject(const CompactSplineIrregular1D& obj, char* newFlatBufferPtr)
{
  /// See FlatObject for description
  const char* oldFlatBufferPtr = obj.mFlatBufferPtr;
  FlatObject::cloneFromObject(obj, newFlatBufferPtr);
  mNumberOfKnots = obj.mNumberOfKnots;
  mUmax = obj.mUmax;
  mBin2KnotMap = FlatObject::relocatePointer(oldFlatBufferPtr, mFlatBufferPtr, obj.mBin2KnotMap);
}

void CompactSplineIrregular1D::moveBufferTo(char* newFlatBufferPtr)
{
  /// See FlatObject for description
  const char* oldFlatBufferPtr = mFlatBufferPtr;
  FlatObject::moveBufferTo(newFlatBufferPtr);
  mBin2KnotMap = FlatObject::relocatePointer(oldFlatBufferPtr, mFlatBufferPtr, mBin2KnotMap);
}
#endif

void CompactSplineIrregular1D::setActualBufferAddress(char* actualFlatBufferPtr)
{
  /// See FlatObject for description
  mBin2KnotMap = FlatObject::relocatePointer(mFlatBufferPtr, actualFlatBufferPtr, mBin2KnotMap);
  FlatObject::setActualBufferAddress(actualFlatBufferPtr);
}

void CompactSplineIrregular1D::setFutureBufferAddress(char* futureFlatBufferPtr)
{
  /// See FlatObject for description
  mBin2KnotMap = FlatObject::relocatePointer(mFlatBufferPtr, futureFlatBufferPtr, mBin2KnotMap);
  FlatObject::setFutureBufferAddress(futureFlatBufferPtr);
}

#if !defined(GPUCA_GPUCODE)

void CompactSplineIrregular1D::construct(int numberOfKnots, const int inputKnots[])
{
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

  FlatObject::startConstruction();

  std::vector<int> knotU;

  { // reorganize knots

    std::vector<int> tmp;
    for (int i = 0; i < numberOfKnots; i++)
      tmp.push_back(inputKnots[i]);
    std::sort(tmp.begin(), tmp.end());

    knotU.push_back(0); // obligatory knot at 0.0

    for (int i = 0; i < numberOfKnots; ++i) {
      if (knotU.back() < tmp[i])
        knotU.push_back(tmp[i]);
    }
    if (knotU.back() < 1)
      knotU.push_back(1);
  }

  mNumberOfKnots = knotU.size();
  mUmax = knotU.back();
  int bin2KnotMapOffset = mNumberOfKnots * sizeof(CompactSplineIrregular1D::Knot);

  FlatObject::finishConstruction(bin2KnotMapOffset + (mUmax + 1) * sizeof(int));

  mBin2KnotMap = reinterpret_cast<int*>(mFlatBufferPtr + bin2KnotMapOffset);

  CompactSplineIrregular1D::Knot* s = getKnotsNonConst();

  for (int i = 0; i < mNumberOfKnots; i++) {
    s[i].u = knotU[i];
  }

  for (int i = 0; i < mNumberOfKnots - 1; i++) {
    s[i].Li = 1. / (s[i + 1].u - s[i].u); // do division in double
  }

  s[mNumberOfKnots - 1].Li = 0.f; // the value will not be used, we define it for consistency

  // Set up map (U bin) -> (knot index)

  int* map = getBin2KnotMapNonConst();

  int iKnotMax = mNumberOfKnots - 2;

  //
  // With iKnotMax=nKnots-2 we map the U==Umax coordinate to the [nKnots-2, nKnots-1] segment.
  // This trick allows one to avoid a special condition for this edge case.
  // Any U from [0,Umax] is mapped to some knot i such, that the knot i+1 is always exist
  //

  for (int u = 0, iKnot = 0; u <= mUmax; u++) {
    if ((knotU[iKnot + 1] == u) && (iKnot < iKnotMax)) {
      iKnot = iKnot + 1;
    }
    map[u] = iKnot;
  }
}

void CompactSplineIrregular1D::constructRegular(int numberOfKnots)
{
  /// Constructor for a regular spline
  /// \param numberOfKnots     Number of knots
  ///

  if (numberOfKnots < 2)
    numberOfKnots = 2;

  std::vector<int> knots(numberOfKnots);
  for (int i = 0; i < numberOfKnots; i++) {
    knots[i] = i;
  }
  construct(numberOfKnots, knots.data());
}
#endif

void CompactSplineIrregular1D::print() const
{
#if !defined(GPUCA_GPUCODE)
  std::cout << " Compact Spline 1D: " << std::endl;
  std::cout << "  mNumberOfKnots = " << mNumberOfKnots << std::endl;
  std::cout << "  mUmax = " << mUmax << std::endl;
  std::cout << "  mBin2KnotMap = " << (void*)mBin2KnotMap << std::endl;
  std::cout << "  knots: ";
  for (int i = 0; i < mNumberOfKnots; i++) {
    std::cout << getKnot(i).u << " ";
  }
  std::cout << std::endl;
#endif
}
