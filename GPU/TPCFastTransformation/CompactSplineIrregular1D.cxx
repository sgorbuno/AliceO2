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

CompactSplineIrregular1D::CompactSplineIrregular1D() : FlatObject(), mNumberOfKnots(0), mNumberOfAxisBins(0), mBin2KnotMapOffset(0)
{
  /// Default constructor. Creates an empty uninitialised object
}

void CompactSplineIrregular1D::destroy()
{
  /// See FlatObject for description
  mNumberOfKnots = 0;
  mNumberOfAxisBins = 0;
  mBin2KnotMapOffset = 0;
  FlatObject::destroy();
}

void CompactSplineIrregular1D::cloneFromObject(const CompactSplineIrregular1D& obj, char* newFlatBufferPtr)
{
  /// See FlatObject for description
  FlatObject::cloneFromObject(obj, newFlatBufferPtr);
  mNumberOfKnots = obj.mNumberOfKnots;
  mNumberOfAxisBins = obj.mNumberOfAxisBins;
  mBin2KnotMapOffset = obj.mBin2KnotMapOffset;
}

void CompactSplineIrregular1D::construct(int numberOfKnots, const float inputKnots[], int numberOfAxisBins)
{
  /// Constructor.
  /// Initializes the spline with a grid with numberOfKnots knots in the interval [0,1)
  /// array inputKnots[] has numberOfKnots entries in increase order, started from 0.
  /// A knot on the edges u==0. is obligatory
  ///
  /// The number of knots created and their values may change during initialisation:
  /// - Edge knot 0.f will be added when it is not present.
  /// - Knot values are rounded to closest axis bins: k*1./numberOfAxisBins.
  /// - Knots which are rounded to the same bin will be merged
  /// - At least 2 knots and at least 2 axis bins will be created
  ///
  /// \param numberOfKnots     Number of knots in knots[] array
  /// \param knots             Array of knots.
  /// \param numberOfAxisBins Number of axis bins to map U coordinate to
  ///                          an appropriate [knot(i),knot(i+1)] interval.
  ///                          The knot positions will have a "granularity" of 1./numberOfAxisBins
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
  mBin2KnotMapOffset = mNumberOfKnots * sizeof(CompactSplineIrregular1D::Knot);

  FlatObject::finishConstruction(mBin2KnotMapOffset + (numberOfAxisBins + 1) * sizeof(int));

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
  // With iKnotMax=nKnots-2 we map U coordinates from the outside of the intreval [nKnots-1, inf) to [nKnots-2, nKnots-1)
  // This trick allows one to use splines without any special condition for the edge case.
  // Any U from [0,1] is mapped to some knot i such, that the knot i+1 is always exist
  //

  for (int iBin = 0, iKnot = 0; iBin <= mNumberOfAxisBins; iBin++) {
    if ((iKnot < iKnotMax) && vKnotBins[iKnot + 1] == iBin) {
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
  std::cout << "  mBin2KnotMapOffset = " << mBin2KnotMapOffset << std::endl;
  std::cout << "  knots: ";
  for (int i = 0; i < mNumberOfKnots; i++) {
    std::cout << getKnot(i).u << " ";
  }
  std::cout << std::endl;
#endif
}
