// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  SplineXD.cxx
/// \brief Implementation of SplineXD class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#include "SplineXD.h"

#include <iostream>

templateClassImp(GPUCA_NAMESPACE::gpu::SplineXD);

using namespace std;
using namespace GPUCA_NAMESPACE::gpu;

template class GPUCA_NAMESPACE::gpu::SplineXD<float>;
template class GPUCA_NAMESPACE::gpu::SplineXD<double>;

//template class GPUCA_NAMESPACE::gpu::SplineXD<float, 0>;
//template class GPUCA_NAMESPACE::gpu::SplineXD<double, 0>;
