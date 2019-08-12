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
 * test_mesh.cpp
 *
 *  Created on: Nov 24, 2014
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
#include <Tpetra_Map.hpp>

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include <string>
#include <vector>
using namespace std;

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
#if 1
TEUCHOS_UNIT_TEST( Int, Mesh_RefToPhysInternal ) {
  Mesh mesh(mesh_file, comm);
  LO ielement = 11;
  VECTOR ref, phys, ref2;
  ref[0] = 0.1; ref[1] = 0.2; ref[2] = 0.3;
  mesh.refToPhysical(phys, ref, ielement);
  mesh.physToReference(ref2, phys, ielement);
  FLOAT delta = 1;
  for (int i=0;i<3;++i)
    delta += (ref[i]-ref2[i])*(ref[i]-ref2[i]);
  TEST_FLOATING_EQUALITY(delta, 1., 1e-12);

}
TEUCHOS_UNIT_TEST( Int, Mesh_RefToPhysExternal ) {
  Mesh mesh(mesh_file, comm);
  LO ielement = 11;
  VECTOR ref, phys, ref2;
  ref[0] = 2.1; ref[1] = -1.2; ref[2] = 3.3;
  mesh.refToPhysical(phys, ref, ielement);
  mesh.physToReference(ref2, phys, ielement);
  FLOAT delta = 1;
  for (int i=0;i<3;++i)
    delta += (ref[i]-ref2[i])*(ref[i]-ref2[i]);
  TEST_FLOATING_EQUALITY(delta, 1., 1e-12);

}


TEUCHOS_UNIT_TEST( Int, Mesh_RefToPhysVectorInternal ) {
  Mesh mesh(mesh_file, comm);
  LO ielement = 11;
  VECTOR ref;
  ref[0] = 0.1; ref[1] = 0.2; ref[2] = 0.3;
//  ref[0] = ref[1] = ref[2] = -0.0;
  VECTOR v, v_ref, v_test;
  v[0] = v[1] = v[2] = 1.0;
  mesh.physToReferenceVector(v_ref, v, ref, ielement);
  mesh.refToPhysicalVector(v_test, v_ref, ref, ielement);
  FLOAT delta = 1;
  for (int i=0;i<3;++i)
    delta += (v[i]-v_test[i])*(v[i]-v_test[i]);
  TEST_FLOATING_EQUALITY(delta, 1., 1e-12);

}


TEUCHOS_UNIT_TEST( Int, Mesh_is_in ) {
  Mesh mesh(mesh_file, comm);
  LO ielement = 11;
  VECTOR ref, phys, ref2;
  ref[0] = 0.1; ref[1] = 0.2; ref[2] = 0.3;
  mesh.refToPhysical(phys, ref, ielement);
  bool is_in = mesh.isPhysicalInElement(phys, ielement);
  TEST_ASSERT(is_in);

  ref[0] = 2.1; ref[1] = -1.2; ref[2] = 3.3;
  mesh.refToPhysical(phys, ref, ielement);
  is_in = mesh.isPhysicalInElement(phys, ielement);
  TEST_ASSERT(!is_in);

}

// Make sure that Mesh::Face works as a key with std::map
TEUCHOS_UNIT_TEST( Int, Face_std_map ) {
  const int DIM = 3;
  typedef Mesh::Face<DIM> Face;

  const int size = 5;
  std::map<Face, int> m;
  std::vector<Face> faces(size);
  for(int i=0; i<size; ++i) {
    ConstLenVector<GO, DIM> v;
    for(int d=0; d<DIM; ++d) v[d] = i;
    faces[i] = Face(v);
    m[faces[i]] = i;
  }

  ConstLenVector<GO, DIM> vs;
  for(int d=0; d<DIM; ++d) vs[d] = 1;
  Face s(vs);

  auto it = m.find(s);
  TEST_ASSERT(it!=m.end());

  m[s] = 999;
  TEST_ASSERT(m.size()==size);
  TEST_ASSERT(m[faces[1]]==999);

}

#endif
TEUCHOS_UNIT_TEST( Int, Mesh_build_periodic ) {
  Mesh mesh(mesh_file, comm, true);
}

TEUCHOS_UNIT_TEST( Int, Mesh_volume ) {
  Mesh mesh(mesh_file, comm);
  FLOAT sum = 0, global_sum;
  for (LO i=0; i<mesh.num_elems; ++i)
    sum += mesh.determinate_jacobian(i);

  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &sum, &global_sum);
  TEST_FLOATING_EQUALITY(global_sum, 48., 1e-6); // Tets have a base area of 1/6

}
