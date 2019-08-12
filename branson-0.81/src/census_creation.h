//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   census_creation.h
 * \author Alex Long
 * \date   January 1 2015
 * \brief  Function for creating initial census particles
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef census_creation_h_
#define census_creation_h_

#include <vector>

#include "photon.h"

double get_photon_list_E(const std::vector<Photon> &photons) {
  double total_E = 0.0;
  for (auto const &iphtn : photons)
    total_E += iphtn.get_E();
  return total_E;
}

#endif // def census_creation_h_
//---------------------------------------------------------------------------//
// end of census_creation.h
//---------------------------------------------------------------------------//
