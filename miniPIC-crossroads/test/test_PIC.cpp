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
 * test_two_particle.cpp
 *
 *  Created on: Nov 20, 2014
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
#include "PIC.hpp"
#include "test_utils.hpp"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include <string>

Teuchos::RCP<Teuchos::Comm<int> >  comm;
std::string mesh_file;
int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  comm = Teuchos::RCP<Teuchos::Comm<int> > (new Teuchos::MpiComm<int>(MPI_COMM_WORLD) );

  Kokkos::initialize(argc, argv);

  mesh_file = std::string("regular.txt");
  int result =  Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  mesh_file = std::string("mapped.txt");
  result +=  Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
  Kokkos::finalize();
  return result;
}

FLOAT soln(VECTOR &phys) {
  return phys[0]*phys[0] + 3*phys[1] + phys[2]*phys[2]*phys[2];
}
VECTOR field(VECTOR &phys) {
  VECTOR E;
  E[0] = 2*phys[0]; E[1] = 3; E[2] = 3*phys[2]*phys[2];
  return E;
}

TEUCHOS_UNIT_TEST( Int, OneParticleAnalyticForce ) {

  DataWarehouse data;
  data.create_mesh(mesh_file, comm);
  data.create_particle_list(1000);
  data.create_ES();
  PIC pic(data);

  data.ES_reference().set_analytic_phi(soln);

  Particle p;
  p.x[0] = -0.11234; p.x[1] = -0.1234; p.x[2] = -0.134;
  p.v[0] = p.v[1] = p.v[2] = 0;
  LO ielement;
  Kokkos::parallel_reduce(data.mesh_reference().num_elems,  FindElement(data.mesh_reference(), p.x), ielement);
  p.ielement = ielement;
  if (p.ielement >= 0 ) {
    data.mesh_reference().physToReference(p.x_ref, p.x, p.ielement);
    FillOneParticle f1(data, p);
    f1.execute();
  }
  pic.weight_Efield();

  if ( data.particles().num_particles()) {
    p.E = field(p.x);
    CompareOnDevice comp(p, data.particles(), 0, 2e-1, CompareOnDevice::E);
    bool result = comp.compare();
    TEST_EQUALITY_CONST(true, result);

  }
}

TEUCHOS_UNIT_TEST( Int, TwoParticlesComputeForce ) {

  DataWarehouse data;
  data.create_mesh(mesh_file, comm);
  data.create_particle_list(1000);
  data.create_ES();
  PIC pic(data);

  Particle p;
  p.type = 0;
  p.x[0] = 0.11234; p.x[1] = 0.1234; p.x[2] = 0.134;
  p.v[0] = p.v[1] = p.v[2] = 0;
  p.weight=1;
  LO ielement;
  Kokkos::parallel_reduce(data.mesh_reference().num_elems,  FindElement(data.mesh_reference(), p.x), ielement);
  p.ielement = ielement;
  if (p.ielement >= 0 ) {
    data.mesh_reference().physToReference(p.x_ref, p.x, p.ielement);
    FillOneParticle f1(data, p);
    f1.execute();
  }
  p.type = 1;
  p.x[0] = -p.x[0]; p.x[1] = -p.x[1]; p.x[2] = -p.x[2];
  p.weight=0;
  Kokkos::parallel_reduce(data.mesh_reference().num_elems,  FindElement(data.mesh_reference(), p.x), ielement);
  if (p.ielement >= 0 ) {
    data.mesh_reference().physToReference(p.x_ref, p.x, p.ielement);
    FillOneParticle f2(data, p);
    f2.execute();
  }

  FLOAT dist = 2*sqrt(p.x[0]*p.x[0]+p.x[1]*p.x[1]+p.x[2]*p.x[2]);
  std::cout << "solution should be "<<1./(4*M_PI*dist*dist) << std::endl;
  pic.time_step(.1);

  std::cout << "E total " << sqrt(data.particles_reference().E(1)[0]*data.particles_reference().E(1)[0]+data.particles_reference().E(1)[1]*data.particles_reference().E(1)[1]+
      data.particles_reference().E(1)[2]*data.particles_reference().E(1)[2])<<endl;
  std::cout << data.particles_reference().E(0)[0] << " "<< data.particles_reference().E(0)[1] << " "<< data.particles_reference().E(0)[2] << " "
            << data.particles_reference().E(1)[0] << " "<< data.particles_reference().E(1)[1] << " "<< data.particles_reference().E(1)[2] << "\n";
}
