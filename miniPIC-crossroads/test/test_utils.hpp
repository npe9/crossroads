/**********************************************************************************

Mini-PIC 

Copyright (c) 2015, Sandia National Laboratories
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

For questions, comments or contributions contact 
Matt Bettencourt, mbetten@sandia.gov

*******************************************************************************/
/*
 * utils.hpp
 *
 *  Created on: Nov 19, 2014
 *      Author: mbetten
 */

#ifndef TEST_UTILS_HPP_
#define TEST_UTILS_HPP_

#include "types.hpp"
#include "mesh.hpp"
#include "particle_list.hpp"
#include "particle_move.hpp"
#include "particle_fill.hpp"
#include "base_functor.hpp"
#include "data_warehouse.hpp"

#include <algorithm>
//namespace {

class FillOneParticle{
public:
  FillOneParticle(DataWarehouse &data, Particle &p):
    data_(data), p_(p), parts_(data_.particles()){}
  void execute(){
    TeamPolicy policy(1, 1);
    Kokkos::fence();
    Kokkos::parallel_for(policy, *this);
    Kokkos::fence();
  }
  KOKKOS_INLINE_FUNCTION
  void operator () (const member_type & dev) const {
    long is_filled=true;
    parts_.push_back(dev, &p_, &is_filled, 1);
  }
  DataWarehouse &data_;
  Particle p_;
  ParticleList parts_;

};

struct FindElement{
  FindElement(Mesh &mesh, VECTOR &p) :mesh_(mesh), p_(p){ }
  typedef LO value_type;
  KOKKOS_INLINE_FUNCTION
  void join(volatile LO &dst, const volatile LO &src ) const {if ( src > dst ) dst = src;}
  KOKKOS_INLINE_FUNCTION
  void init(LO &v)const {v = -1;}
  KOKKOS_INLINE_FUNCTION
  void operator()(std::size_t indx, LO &update ) const {
    if (mesh_.isPhysicalInElement(p_, indx))
      update = indx;
   }
  Mesh mesh_;
  VECTOR p_;
};

struct CompareOnDevice{
  enum which_field {x=0, v, E};

  CompareOnDevice(Particle p, ParticleList parts, std::size_t index, FLOAT tol, which_field field= x ) :
    p_(p), parts_(parts), index_(index), tol_(tol), field_(field), result_("Result") {}

  bool compare(){
    Kokkos::fence();
    Kokkos::parallel_for(1, *this);
    Kokkos::deep_copy(result_.h_view, result_.d_view);
    Kokkos::fence();

    return result_.h_view();
  }

  KOKKOS_INLINE_FUNCTION
  void operator ()(int i) const {
    result_.d_view() = true;
    if ( parts_.num_particles() > index_)
      for (int i=0; i<3; ++i ) 
        if ( field_ == x )
          result_.d_view() = ( result_.d_view() && p_.x[i] + tol_ > parts_.x(index_)[i] && p_.x[i] - tol_ < parts_.x(index_)[i] );
        else if ( field_ == v )
          result_.d_view() = ( result_.d_view() && p_.v[i] + tol_ > parts_.v(index_)[i] && p_.v[i] - tol_ < parts_.v(index_)[i] );
        else if ( field_ == E )
          result_.d_view() = ( result_.d_view() && p_.E[i] + tol_ > parts_.E(index_)[i] && p_.E[i] - tol_ < parts_.E(index_)[i] );
        else
          result_.d_view() = false;

  }
  Particle p_;
  ParticleList parts_;
  size_t index_;
  FLOAT tol_;
  which_field field_;
  Kokkos::DualView<bool> result_;
};


inline
void CreateOneParticle(DataWarehouse & data, std::string &file, Particle &p, Teuchos::RCP<Teuchos::Comm<int> >  comm, bool periodic=false) {
   data.create_mesh(file, comm, periodic);
   data.create_particle_list(1234);

   p.ielement = 1; // This is GID of 2
   p.x_ref[0] = p.x_ref[1] = p.x_ref[2] = .1234;
   data.mesh().refToPhysical(p.x, p.x_ref, p.ielement);

   p.v[0] = .5; p.v[1] = .25; p.v[2] = .15;
   p.is_dead = false;

   if (data.mesh().element_map->isNodeGlobalElement(p.ielement+1)) {
     FillOneParticle fill(data, p);
     fill.execute();
   }
}



#endif /* TEST_UTILS_HPP_ */
