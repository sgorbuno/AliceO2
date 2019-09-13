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
#include <vector>

#if !defined(GPUCA_GPUCODE)
#include <iostream>
#endif

using namespace GPUCA_NAMESPACE::gpu;

CompactSplineIrregular1D::CompactSplineIrregular1D() : FlatObject(), mNumberOfKnots(0), mNumberOfAxisBins(0), mBin2KnotMap(0)
{
  /// Default constructor. Creates an empty uninitialised object
}

void CompactSplineIrregular1D::destroy()
{
  /// See FlatObject for description
  mNumberOfKnots = 0;
  mNumberOfAxisBins = 0;
  mBin2KnotMap = nullptr;
  FlatObject::destroy();
}

void CompactSplineIrregular1D::cloneFromObject(const CompactSplineIrregular1D& obj, char* newFlatBufferPtr)
{
  /// See FlatObject for description
  FlatObject::cloneFromObject(obj, newFlatBufferPtr);
  mNumberOfKnots = obj.mNumberOfKnots;
  mNumberOfAxisBins = obj.mNumberOfAxisBins;
  mBin2KnotMap = FlatObject::relocatePointer(obj.mFlatBufferPtr, mFlatBufferPtr, obj.mBin2KnotMap);
}

void IrregularSpline2D3D::moveBufferTo(char* newFlatBufferPtr)
{
/// See FlatObject for description
#ifndef GPUCA_GPUCODE
  FlatObject::moveBufferTo(newFlatBufferPtr);
  setActualBufferAddress(mFlatBufferPtr);
#endif
}

void IrregularSpline2D3D::setActualBufferAddress(char* actualFlatBufferPtr)
{
  /// See FlatObject for description
  mBin2KnotMap = FlatObject::relocatePointer(mFlatBufferPtr, actualFlatBufferPtr, mBin2KnotMap);
  FlatObject::setActualBufferAddress(actualFlatBufferPtr);
}

void IrregularSpline2D3D::setFutureBufferAddress(char* futureFlatBufferPtr)
{
  /// See FlatObject for description
  mBin2KnotMap = FlatObject::relocatePointer(mFlatBufferPtr, futureFlatBufferPtr, mBin2KnotMap);
  FlatObject::setFutureBufferAddress(futureFlatBufferPtr);
}

void CompactSplineIrregular1D::construct(int numberOfKnots, const float inputKnots[], int numberOfAxisBins)
{
  /// Constructor.
  /// Initializes the spline with a grid with numberOfKnots knots in the segment [0,1)
  /// array inputKnots[] has numberOfKnots entries in increase order, started from 0.
  /// The edge knots u==0. and u==1. are obligatory
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

  FlatObject::startConstruction();

  if (numberOfAxisBins < 1) {
    numberOfAxisBins = 1;
  }

  std::vector<int> vKnotBins;

  { // reorganize knots

    int lastBin = numberOfAxisBins; // last bin starts with U value 1.f, therefore it is outside of the [0.,1.] interval

    vKnotBins.push_back(0); // obligatory knot at 0.0

    for (int i = 0; i < numberOfKnots; ++i) {
      int bin = (int)roundf(inputKnots[i] * numberOfAxisBins);
      if (bin <= vKnotBins.back() || bin >= lastBin) {
        continue; // same knot
      }
      vKnotBins.push_back(bin);
    }

    vKnotBins.push_back(lastBin); // obligatory knot at 1.0
  }

  mNumberOfKnots = vKnotBins.size();
  mNumberOfAxisBins = numberOfAxisBins;
  int bin2KnotMapOffset = mNumberOfKnots * sizeof(CompactSplineIrregular1D::Knot);

  FlatObject::finishConstruction(bin2KnotMapOffset + (numberOfAxisBins + 1) * sizeof(int));

  mBin2KnotMap = reinterpret_cast<int*>(mFlatBufferPtr + bin2KnotMapOffset);

  CompactSplineIrregular1D::Knot* s = getKnotsNonConst();

  for (int i = 0; i < mNumberOfKnots; i++) {
    s[i].u = vKnotBins[i] / ((double)mNumberOfAxisBins); // do division in double
  }

  for (int i = 0; i < mNumberOfKnots - 1; i++) {
    s[i].L = s[i + 1].u - s[i].u;
    s[i].Li = 1. / s[i].L;
  }

  { // values will not be used, we define them for consistency
    int i = mNumberOfKnots - 1;
    s[i].L = 0.f;
    s[i].Li = 0.f;
  }

  // Set up map (U bin) -> (knot index)

  int* map = getBin2KnotMapNonConst();

  int iKnotMax = mNumberOfKnots - 2;

  //
  // With iKnotMax=nKnots-2 we map the U==1 coordinate to the [nKnots-2, nKnots-1] segment.
  // This trick allows one to avoid a special condition for the edge case.
  // Any U from [0,1] is mapped to some knot i such, that the knot i+1 is always exist
  //

  for (int iBin = 0, iKnot = 0; iBin <= mNumberOfAxisBins; iBin++) {
    if ((vKnotBins[iKnot + 1] == iBin) && (iKnot < iKnotMax)) {
      iKnot = iKnot + 1;
    }
    map[iBin] = iKnot;
  }
}

void CompactSplineIrregular1D::constructRegular(int numberOfKnots)
{
  /// Constructor for a regular spline
  /// \param numberOfKnots     Number of knots
  ///

  if (numberOfKnots < 2)
    numberOfKnots = 2;

  std::vector<float> knots(numberOfKnots);
  double du = 1. / (numberOfKnots - 1);
  for (int i = 0; i < numberOfKnots; i++) {
    knots[i] = i * du;
  }
  construct(numberOfKnots, knots.data(), numberOfKnots);
}

void CompactSplineIrregular1D::print() const
{
#if !defined(GPUCA_GPUCODE)
  std::cout << " Compact Spline 1D: " << std::endl;
  std::cout << "  mNumberOfKnots = " << mNumberOfKnots << std::endl;
  std::cout << "  mNumberOfAxisBins = " << mNumberOfAxisBins << std::endl;
  std::cout << "  mBin2KnotMap = " << (void*)mBin2KnotMap << std::endl;
  std::cout << "  knots: ";
  for (int i = 0; i < mNumberOfKnots; i++) {
    std::cout << getKnot(i).u << " ";
  }
  std::cout << std::endl;
#endif
}
