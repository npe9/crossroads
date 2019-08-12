//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   photon.h
 * \author Alex Long
 * \date   July 18 2014
 * \brief  Holds values and functions needed for transporting photon
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//----------------------------------------------------------------------------//

#ifndef photon_h_
#define photon_h_

#include <cmath>
#include <iostream>
#include <vector>

#include "constants.h"

//==============================================================================
/*!
 * \class Photon
 * \brief Contains position, direction, cell ID and energy for transport.
 *
 * Holds all of the internal state of a photon and provides functions for
 * sorting photons based on census and global cell ID.
 */
//==============================================================================
class Photon {
public:
  //! Constructor
  Photon() {}

  //! Destructor
  ~Photon(void) {}

  //--------------------------------------------------------------------------//
  // const functions                                                          //
  //--------------------------------------------------------------------------//

  //! Check to see if photon energy weight is below cutoff fraction
  bool below_cutoff(const double cutoff_fraction) const {
    return (m_E / m_E0 < cutoff_fraction);
  }

  //! Return global cell ID
  inline uint32_t get_cell(void) const { return m_cell_ID; }

  //! Return global grip ID
  inline uint32_t get_grip(void) const { return m_grip_ID; }

  //! Return photon group
  inline uint32_t get_group(void) const { return group; }

  //! Return a constant pointer to the start of the particle position array
  inline const double *get_position(void) const { return m_pos; }

  //! Return a constant pointer to the start of the particle direction array
  inline const double *get_angle(void) const { return m_angle; }

  //! Get the particle's energy-weight
  inline double get_E(void) const { return m_E; }

  //! Get the partice's initial energy-weight
  inline double get_E0(void) const { return m_E0; }

  //! Get the distance to census (cm)
  inline double get_distance_remaining(void) const { return m_life_dx; }

  //! Print particle information
  void print_info(const uint32_t &rank) const {
    using std::cout;
    using std::endl;
    cout << "----Photon Info----\n";
    cout << rank << " " << m_pos[0] << " " << m_pos[1] << " " << m_pos[2]
         << endl;
    cout << "angle: " << m_angle[0] << " " << m_angle[1] << " " << m_angle[2]
         << endl;
    cout << "Energy: " << m_E << " , Initial energy: " << m_E0 << endl;
    cout << "Cell ID: " << m_cell_ID << " , Grip ID: " << m_grip_ID << endl;
  }

  //! Override great than operator to sort
  bool operator<(const Photon &compare) const {
    return m_cell_ID < compare.get_cell();
  }

  //--------------------------------------------------------------------------//
  // non-const functions                                                      //
  //--------------------------------------------------------------------------//

  //! Update particle position by moving it some distance
  inline void move(const double distance) {
    m_pos[0] += m_angle[0] * distance;
    m_pos[1] += m_angle[1] * distance;
    m_pos[2] += m_angle[2] * distance;
    m_life_dx -= distance;
  }

  //! Set the global cell ID
  inline void set_cell(const uint32_t new_cell) { m_cell_ID = new_cell; }

  //! Set the global grip ID
  inline void set_grip(const uint32_t new_grip) { m_grip_ID = new_grip; }

  //! Set the group of the photon
  inline void set_group(const uint32_t new_group) { group = new_group; }

  //! Set the initial energy-weight
  inline void set_E0(const double E) {
    m_E0 = E;
    m_E = E;
  }

  //! Set the current energy-weight
  inline void set_E(const double E) { m_E = E; }

  //! Set the distance to census (cm)
  inline void set_distance_to_census(const double dist_remain) {
    m_life_dx = dist_remain;
  }

  //! Set the angle of the photon
  inline void set_angle(double const *const angle) {
    m_angle[0] = angle[0];
    m_angle[1] = angle[1];
    m_angle[2] = angle[2];
  }

  //! Set the spatial position of the photon
  inline void set_position(double const *const pos) {
    m_pos[0] = pos[0];
    m_pos[1] = pos[1];
    m_pos[2] = pos[2];
  }

  //! Reflect a photon about a plane aligned with the X, Y, or Z axes
  inline void reflect(const uint32_t surface_cross) {
    using Constants::X_NEG;
    using Constants::X_POS;
    using Constants::Y_NEG;
    using Constants::Y_POS;
    //reflect the photon over the surface it was crossing
    if (surface_cross == X_POS || surface_cross == X_NEG)
      m_angle[0] = -m_angle[0];
    else if (surface_cross == Y_POS || surface_cross == Y_NEG)
      m_angle[1] = -m_angle[1];
    else
      m_angle[2] = -m_angle[2];
  }

  //--------------------------------------------------------------------------//
  // member data                                                              //
  //--------------------------------------------------------------------------//
private:
  uint32_t m_cell_ID; //!< Cell ID
  uint32_t m_grip_ID; //!< Grip ID of current cell
  uint32_t group;     //!< Group of photon
  double m_pos[3];    //!< photon position
  double m_angle[3];  //!< photon angle array
  double m_E;         //!< current photon energy
  double m_E0;        //!< photon energy at creation
  double m_life_dx;   //!< Distance remaining this time step

  // private member functions
private:
};

#endif // photon_h_
