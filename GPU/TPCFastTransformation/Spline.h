// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  Spline.h
/// \brief Definition of Spline class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#ifndef ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SPLINE_H
#define ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_SPLINE_H

#include "SplineSpec.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
///
/// The Spline class performs a cubic spline interpolation on an two-dimensional nonunifom grid.
/// The class is an extension of the Spline1D class.
/// See Spline1D.h for more details.
///
/// The spline S(x1,x2) approximates a function F(x1,x2):R^2->R^m,
/// with 2-dimensional domain and multi-dimensional codomain.
/// x1,x2 belong to [x1min,x1max] x [x2min,x2max].
///
/// --- Example of creating a spline ---
///
///  auto F = [&](float x1, float x2, float f[] ) {
///   f[0] = 1.f + x1 + x2*x2; // F(x1,x2)
///  };
///  const int nKnotsU=2;
///  const int nKnotsV=3;
///  int knotsU[nKnotsU] = {0, 1};
///  int knotsV[nKnotsV] = {0, 2, 5};
///  Spline<float,1> spline(nKnotsU, knotsU, nKnotsV, knotsV ); // spline with 1-dimensional codomain
///  spline.approximateFunction(0., 1., 0.,1., F); //initialize spline to approximate F on [0., 1.]x[0., 1.] area
///  float S = spline.interpolate(.1, .3 ); // interpolated value at (.1,.3)
///
///  --- See also SplineHelper::test();
///

/// ==================================================================================================
///
/// Declare the Spline class as a template with one optional parameter.
///
/// The default value is just an indicator of the absence of the second parameter.
/// (The right way would be to use variadic templates for this case,
/// but they are screwed up in the ROOT linker).
///
/// Class specifications depend on the YdimT value. They can be found in SplineSpecs.h
///
/// \param DataT data type: float or double
/// \param YdimT >= 0 : number of Y dimensions,
///               < 0 : max possible number of Y dimensions
///              default : no info about Y dimensions
///
template <typename DataT, int XdimT = 0, int YdimT = 0>
class Spline
  : public SplineSpec<DataT,
                      XdimT, (XdimT > 0), false,
                      YdimT, (YdimT > 0), false>
{
  typedef SplineContainer<DataT> TVeryBase;
  typedef SplineSpec<DataT,
                     XdimT, (XdimT > 0), false,
                     YdimT, (YdimT > 0), false>
    TBase;

 public:
  typedef typename TVeryBase::SafetyLevel SafetyLevel;
  typedef typename TVeryBase::Knot Knot;

#if !defined(GPUCA_GPUCODE)
  using TBase::TBase; // inherit constructors

  /// Assignment operator
  Spline& operator=(const Spline& v)
  {
    TVeryBase::cloneFromObject(v, nullptr);
    return *this;
  }
#else
  /// Disable constructors for the GPU implementation
  Spline() CON_DELETE;
  Spline(const SplineSpec&) CON_DELETE;
#endif

  /// cast to the one-parameter-template class
  operator Spline<DataT>&() { return *((Spline<DataT>*)this); }

  /// const cast
  operator const Spline<DataT>&() const { return *((const Spline<DataT>*)this); }

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE)
  /// read a class object from the file
  static Spline* readFromFile(TFile& inpf, const char* name)
  {
    return (Spline*)TVeryBase::readFromFile(inpf, name);
  }
#endif

#ifndef GPUCA_ALIROOT_LIB
  ClassDefNV(Spline, 0);
#endif
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
