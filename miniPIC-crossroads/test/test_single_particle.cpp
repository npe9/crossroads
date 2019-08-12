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
 * test_single_particle.cpp
 *
 *  Created on: Oct 6, 2014
 *      Author: mbetten
 */




#include <mpi.h>
#include "Kokkos_Core.hpp"
#include <impl/Kokkos_Timer.hpp>
#include "Teuchos_GlobalMPISession.hpp"
#include "types.hpp"
#include "mesh.hpp"
#include "particle_list.hpp"
#include "particle_move.hpp"
#include "particle_fill.hpp"
#include "base_functor.hpp"
#include "data_warehouse.hpp"
#include "test_utils.hpp"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include <string>

Teuchos::RCP<Teuchos::Comm<int> >  comm;
std::string mesh_file;
int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  comm = Teuchos::RCP<Teuchos::Comm<int> > (new Teuchos::MpiComm<int>(MPI_COMM_WORLD) );
  int result;
  {
    Kokkos::initialize(argc, argv);
    mesh_file = std::string("regular.txt");
    result =  Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
    Kokkos::finalize();
  }

  {
    Kokkos::initialize(argc, argv);
    mesh_file = std::string("mapped.txt");
    result =  Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
    Kokkos::finalize();
  }

  return result;
}

TEUCHOS_UNIT_TEST( Int, SingleParticle_SingleStep ) {
  DataWarehouse data;
  Particle p;
  CreateOneParticle(data, mesh_file, p, comm);

  VECTOR initial = p.x;
  VECTOR v = p.v;

  ParticleMove mover(data);
  FLOAT dt = 0.1;
  mover.set_force_scale(0.0);
  Kokkos::fence();
  mover.move(dt);
  Kokkos::fence();

  for (int i=0;i<3; ++i) 
    p.x[i] = initial[i]+dt*v[i];

  if ( data.particles().num_particles()) {
    CompareOnDevice comp(p, data.particles(), 0, 1e-12);
    bool result = comp.compare();
    TEST_EQUALITY_CONST(true, result);
  }
}

#if 1
TEUCHOS_UNIT_TEST( Int, SingleParticle_MultiStep ) {
  DataWarehouse data;
  Particle p;
  CreateOneParticle(data, mesh_file, p, comm);

  VECTOR initial = p.x;
  VECTOR v = p.v;

  ParticleMove mover(data);
  FLOAT dt = 0.1;
  mover.set_force_scale(0.0);
  int nsteps = 10;
  for (int i=0; i<nsteps; ++i)
    mover.move(dt);

  for (int i=0;i<3; ++i) 
    p.x[i] = initial[i]+dt*nsteps*v[i];

  if ( data.particles().num_particles()) {
    CompareOnDevice comp(p, data.particles(), 0, 1e-12);
    bool result = comp.compare();
    TEST_EQUALITY_CONST(true, result);
  }
}


TEUCHOS_UNIT_TEST( Int, SingleParticle_MultiStep_Reflection ) {
  DataWarehouse data;
  Particle p;
  CreateOneParticle(data, mesh_file, p, comm);
  VECTOR initial = p.x;
  VECTOR v = p.v;

  ParticleMove mover(data);
  FLOAT dt = 0.4;
  mover.set_force_scale(0.0);
  int nsteps = 10;
  for (int i=0; i<nsteps; ++i)
    mover.move(dt);

  for (int idim=0;idim<3; ++idim) {
    p.x[idim] = initial[idim]+dt*nsteps*v[idim];
    if ( p.x[idim] > 1)
      p.x[idim] = 1 - (p.x[idim] - 1);
  }
  if ( data.particles().num_particles()) {
    CompareOnDevice comp(p, data.particles(), 0, 1e-12);
    bool result = comp.compare();
    TEST_EQUALITY_CONST(true, result);
  }
}


TEUCHOS_UNIT_TEST( Int, SingleParticle_MultiStep_Periodic ) {
  DataWarehouse data;
  Particle p;
  CreateOneParticle(data, mesh_file, p, comm, true);
  VECTOR initial = p.x;
  VECTOR v = p.v;

  ParticleMove mover(data);
  FLOAT dt = 0.4;
  mover.set_force_scale(0.0);
  int nsteps = 10;
  for (int i=0; i<nsteps; ++i)
    mover.move(dt);

  for (int idim=0;idim<3; ++idim) {
    p.x[idim] = initial[idim]+dt*nsteps*v[idim];
    if ( p.x[idim] > 1)
      p.x[idim] -= 2;
  }
  if ( data.particles().num_particles()) {
    CompareOnDevice comp(p, data.particles(), 0, 1e-12);
    bool result = comp.compare();
    TEST_EQUALITY_CONST(true, result);
  }
}
//}
#endif
