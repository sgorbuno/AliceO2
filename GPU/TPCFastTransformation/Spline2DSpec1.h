// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  Spline2DSpec.h
/// \brief Definition of Spline2DSpec class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#ifndef MYTEST_H
#define MYTEST_H

#include "TObject.h"

template <typename DataTT>
class MyTest : public TObject //FlatObject
{
 public:
  int mYdim;

  ClassDefNV(MyTest, 0);
};

#endif
