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
#include "ES.hpp"
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
  return sin(M_PI*phys[0])*cos(0.5*M_PI*phys[1])*cos(0.5*M_PI*phys[2]);
}
FLOAT rhs(VECTOR &phys) {
  return -1.5*M_PI*M_PI*soln(phys);
}

TEUCHOS_UNIT_TEST( Int, TestAnalyticEField ) {
  Mesh mesh(mesh_file, comm);
  ES es(mesh);
  es.set_analytic_RHS(rhs);
  es.solve();

  FLOAT local_max_error = 0, max_error;
  int num_owned_nodes= mesh.owned_node_map->getNodeNumElements();
  for (unsigned inode=0; inode < num_owned_nodes; ++inode) {
    LO local_node = mesh.node_map->getLocalElement(mesh.owned_node_map->getGlobalElement(inode));
    VECTOR &phys = mesh.nodes(local_node);
    FLOAT p = soln(phys);
    FLOAT s = es.phi(inode);
    local_max_error = std::max(local_max_error, std::abs(s-p));
  }
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &local_max_error, &max_error);
  if ( comm->getRank() == 0 )
    std::cout << "Max error = "<<local_max_error<<endl;
  max_error += 1.0; // normalize the error
  TEST_FLOATING_EQUALITY(1.0, max_error, 0.15);
}


FLOAT shifted_soln(VECTOR &phys) {
  return sin(M_PI*phys[0]+.444)*cos(0.5*M_PI*phys[1])*cos(0.5*M_PI*phys[2]);
}
FLOAT shifted_rhs(VECTOR &phys) {
  return -1.5*M_PI*M_PI*shifted_soln(phys);
}

#if 1
TEUCHOS_UNIT_TEST( Int, TestPeriodicAnalyticEField ) {
  Mesh mesh(mesh_file, comm, true);
  ES es(mesh);
  es.set_analytic_RHS(shifted_rhs);
  es.solve();
  ofstream outs("ES.out", ios::out);
  FLOAT local_max_error = 0, max_error;
  int num_owned_nodes= mesh.owned_node_map->getNodeNumElements();
  for (unsigned inode=0; inode < num_owned_nodes; ++inode) {
    LO local_node = mesh.node_map->getLocalElement(mesh.owned_node_map->getGlobalElement(inode));
    VECTOR &phys = mesh.nodes(local_node);
    FLOAT p = shifted_soln(phys);
    FLOAT s = es.phi(inode);
    local_max_error = std::max(local_max_error, std::abs(s-p));
  }
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &local_max_error, &max_error);
  if ( comm->getRank() == 0 )
    std::cout << "Max error = "<<local_max_error<<endl;
  max_error += 1.0; // normalize the error
  TEST_FLOATING_EQUALITY(1.0, max_error, 0.15);
}
#endif
