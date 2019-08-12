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
#include <Tpetra_Map.hpp>

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include <string>
#include <vector>
using namespace std;

Teuchos::RCP<Teuchos::Comm<int> >  comm;

int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  comm = Teuchos::RCP<Teuchos::Comm<int> > (new Teuchos::MpiComm<int>(MPI_COMM_WORLD) );

  Kokkos::initialize(argc, argv);

  int result =  Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  Kokkos::finalize();
  return result;
}
//namespace {

struct StoreSingleParticle{
  StoreSingleParticle(ParticleList p):parts_(p){}
  KOKKOS_INLINE_FUNCTION
  void operator () (const member_type & dev) const {
    Particle p;
    p.ielement = 123;
    p.is_dead = false;
    long valid = true;
    parts_.push_back(dev, &p, &valid, 1);
  }

  ParticleList parts_;
};

TEUCHOS_UNIT_TEST( Int, StoreSingleParticle ) {

  ParticleList parts(1024);
  
  TeamPolicy policy(1, 1);
  Kokkos::parallel_for(policy, StoreSingleParticle(parts));

  TEST_EQUALITY_CONST(parts.num_particles(), 1);
}


struct StoreMultipleParticle{
  typedef DeviceSpace::scratch_memory_space shmem_space ;
  StoreMultipleParticle(ParticleList p, int nparts, int nelems):parts_(p), nparts_(nparts), nelems_(nelems){}
  unsigned team_shmem_size(int team_size) const {
    return Kokkos::View<Particle* , shmem_space, Kokkos::MemoryUnmanaged>::shmem_size(team_size) +
           Kokkos::View<long*, shmem_space, Kokkos::MemoryUnmanaged>::shmem_size(team_size);
  }
  KOKKOS_INLINE_FUNCTION
  void operator () (const member_type & dev) const {
    const LO leag_rank = dev.league_rank();
    const LO team_size = dev.team_size();
    const LO team_rank = dev.team_rank();
    Kokkos::View<Particle*, shmem_space, Kokkos::MemoryUnmanaged> p(dev.team_shmem(), team_size);
    Kokkos::View<long*, shmem_space, Kokkos::MemoryUnmanaged> is_filled(dev.team_shmem(), team_size);

    GO i = leag_rank*team_size+team_rank;

    p(team_rank).is_dead = false;
    p(team_rank).ielement = ( i % nelems_ );
    p(team_rank).weight = p(team_rank).ielement +100;
    is_filled(team_rank) = ( i < nparts_);
    parts_.push_back(dev, &p(0), &is_filled(0), team_size);
  }

  ParticleList parts_;
  GO nparts_, nelems_;
};
TEUCHOS_UNIT_TEST( Int, StoreMultiParticle1 ) {

  ParticleList parts(1024);

  int nparts = 123;
  StoreMultipleParticle store_parts(parts, nparts, 10);
  int nteams = TeamPolicy::team_size_max(store_parts);
  TeamPolicy policy(nparts/nteams+1, nteams);
  Kokkos::parallel_for(policy, store_parts);

  TEST_EQUALITY_CONST(parts.num_particles(), 123);
}

TEUCHOS_UNIT_TEST( Int, SortParticles ) {

  ParticleList parts(1024);

  int nparts = 123;
  int num_elems = 10;
  StoreMultipleParticle store_parts(parts, nparts, num_elems);
  int nteams = TeamPolicy::team_size_max(store_parts);
  TeamPolicy policy(nparts/nteams+1, nteams);
  Kokkos::parallel_for(policy, store_parts);

  Teuchos::RCP<Map> map = Teuchos::RCP<Map>(new Map(num_elems, 0, comm));
  parts.set_element_map(map);
  bool before_unsort=false, after_unsort=false;
  LO num_parts = parts.num_particles();
  for (LO i=1; i<num_parts; ++i) {
    before_unsort |= (parts.ielement(i-1)> parts.ielement(i));
    before_unsort |= (parts.weight(i) != parts.ielement(i)+100);
  }
  parts.sort();
  for (LO i=1; i<num_parts; ++i){
    after_unsort |=  (parts.ielement(i-1)> parts.ielement(i));
    after_unsort |=  (parts.weight(i) != parts.ielement(i)+100);
  }
  TEST_EQUALITY_CONST(false, after_unsort);
}

TEUCHOS_UNIT_TEST( Int, CountOneParticle ) {

  ParticleList parts(1024);

  int nparts = 1;
  int num_elems = 10;
  StoreMultipleParticle store_parts(parts, nparts, num_elems);
  int nteams = TeamPolicy::team_size_max(store_parts);
  TeamPolicy policy(nparts/nteams+1, nteams);
  Kokkos::parallel_for(policy, store_parts);
  parts.ielement(0) = 4;

  Teuchos::RCP<Map> map = Teuchos::RCP<Map>(new Map(num_elems, 0, comm));
  parts.set_element_map(map);

  parts.compute_start_particles_by_element();
  bool count_right = true;
  for (int i=0; i<num_elems; ++i) {
    LO start = parts.start(i);
    LO count= parts.count_for_element(i);
    if (i <= 3 )
      count_right &= (count == 0 && start == 0 );
    else if (i == 4 )
      count_right &= (count == 1 && start == 0 );
    else
      count_right &= (count == 0 && start == 1 );
  }
  TEST_EQUALITY_CONST(true, count_right);
}

TEUCHOS_UNIT_TEST( Int, CountParticles ) {

  ParticleList parts(1024);

  int nparts = 123;
  int num_elems = 10;
  StoreMultipleParticle store_parts(parts, nparts, num_elems);
  int nteams = TeamPolicy::team_size_max(store_parts);
  TeamPolicy policy(nparts/nteams+1, nteams);
  Kokkos::parallel_for(policy, store_parts);

  Teuchos::RCP<Map> map = Teuchos::RCP<Map>(new Map(num_elems, 0, comm));
  parts.set_element_map(map);

  parts.compute_start_particles_by_element();
  bool start_right = true, count_right = true;
  LO running_start = 0;
  for (int i=0; i<num_elems; ++i) {
    LO start = parts.start(i);
    LO count= parts.count_for_element(i);
    if (i < 3 )
      count_right &= count == 13;
    else
      count_right &= count == 12;
    start_right &= start == running_start;
    running_start += count;
  }
  TEST_EQUALITY_CONST(true, count_right);
  TEST_EQUALITY_CONST(true, start_right);
}


struct MarkForDeletion {
  typedef DeviceSpace::scratch_memory_space shmem_space ;
  MarkForDeletion(ParticleList p, int elem) : parts_(p), which_elem_(elem){}
  unsigned team_shmem_size(int team_size) const {
     return Kokkos::View<std::size_t*, shmem_space, Kokkos::MemoryUnmanaged>::shmem_size(team_size) +
         Kokkos::View<std::size_t, shmem_space, Kokkos::MemoryUnmanaged>::shmem_size();
   }
  KOKKOS_INLINE_FUNCTION
  void operator() (const member_type & dev) const {
   std::size_t i = dev.league_rank() * dev.team_size() + dev.team_rank();
   Kokkos::View<std::size_t*, shmem_space, Kokkos::MemoryUnmanaged> deleted_list(dev.team_shmem(), dev.team_size());
   Kokkos::View<std::size_t, shmem_space, Kokkos::MemoryUnmanaged> count(dev.team_shmem(), dev.team_size());

   if (dev.team_rank() == 0 )
     count() = 0;
   dev.team_barrier();
   if (i < parts_.num_particles() && parts_.ielement(i) == which_elem_ ) {
     std::size_t index = Kokkos::atomic_fetch_add(&count(),1);
     parts_.is_dead(i) = true;
     deleted_list(index) = i;
   }
   dev.team_barrier();
   parts_.delete_particle(dev, &deleted_list(0), count());
 }
  ParticleList parts_;
  GO which_elem_;
};

TEUCHOS_UNIT_TEST( Int, TestDeletion ) {

  ParticleList parts(1024);

  int nparts = 123;
  StoreMultipleParticle store_parts(parts, nparts, 10);
  int nteams = TeamPolicy::team_size_max(store_parts);
  TeamPolicy policy(nparts/nteams+1, nteams);
  Kokkos::parallel_for(policy, store_parts);

  Kokkos::parallel_for(policy, MarkForDeletion(parts, 2));
  parts.cleanup_list();

  TEST_EQUALITY_CONST(parts.num_particles(), 110);

  Kokkos::parallel_for(policy, MarkForDeletion(parts, 5));
  parts.cleanup_list();

  TEST_EQUALITY_CONST(parts.num_particles(), 98);

}
//}

TEUCHOS_UNIT_TEST( Int, TestDeletionBig ) {

  ParticleList parts(1024*128);

  int nparts = 123456;
  StoreMultipleParticle store_parts(parts, nparts, 10);
  int nteams = TeamPolicy::team_size_max(store_parts);
  TeamPolicy policy(nparts/nteams+1, nteams);
  Kokkos::parallel_for(policy, store_parts);

  Kokkos::parallel_for(policy, MarkForDeletion(parts, 2));
  parts.cleanup_list();

  TEST_EQUALITY_CONST(parts.num_particles(), 111110);

  Kokkos::parallel_for(policy, MarkForDeletion(parts, 5));
  parts.cleanup_list();

  TEST_EQUALITY_CONST(parts.num_particles(), 98764);

}

struct FillMigrateList {
  FillMigrateList(Kokkos::View<Kokkos::pair<GO, int>*> migrate_list, int nparts, int nprocs):
    migrate_list_(migrate_list),nparts_(nparts),nprocs_(nprocs), my_rank_(comm->getRank()), comm_size_(comm->getSize()){  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const member_type & dev) const {
    const LO i = dev.league_rank() ;
    LO iproc = i%nprocs_;
    if ( comm_size_ > 1 && comm_size_ == nprocs_+1) {
      // Set a rank other than myown
      if ( iproc >= my_rank_)
        iproc++;
    }
    Kokkos::pair<GO, int> iter((i+91)%nprocs_, iproc);
    migrate_list_(i) = iter;
  }
  Kokkos::View<Kokkos::pair<GO, int>*> migrate_list_;
  int nparts_, nprocs_, my_rank_, comm_size_;
};

TEUCHOS_UNIT_TEST( Int, TestCountList ) {

  int nprocs=13;
  int nparts=123;

  Kokkos::DualView<LO *> proc_reference_list_("ProcList",nprocs);
  Kokkos::View<Kokkos::pair<GO, int>*> migrate_list("MigrateList",1024);

  for (int i=0; i<nprocs; ++i)
    proc_reference_list_.h_view(i) = i;
  Kokkos::deep_copy(proc_reference_list_.d_view,proc_reference_list_.h_view);

  TeamPolicy policy(nparts, 1);
  Kokkos::parallel_for(policy, FillMigrateList(migrate_list, nparts, nprocs));
  Kokkos::fence();

  LO count[nprocs];
  Kokkos::parallel_reduce(nparts,   CountParticleByProcessor(proc_reference_list_.d_view, migrate_list), count);
  Kokkos::fence();
  for (int i=0; i<nprocs; ++i)
    if ( i < 6) {
      TEST_EQUALITY_CONST(count[i], 10);
    }else{
      TEST_EQUALITY_CONST(count[i], 9);
    }

}

PackParticlesForMigration create_packer(int nprocs, int nparts) {
  ParticleList parts(1024);
  StoreMultipleParticle store_parts(parts, nparts, 10);
  int nteams = TeamPolicy::team_size_max(store_parts);
  TeamPolicy policy(nparts/nteams+1, nteams);
  Kokkos::parallel_for(policy, store_parts);


  Kokkos::DualView<LO *> proc_reference_list_("ProcList",nprocs);
  Kokkos::View<Kokkos::pair<GO, int>*> migrate_list("MigrateList",1024);

  int comm_size = comm->getSize();
  int rank = comm->getRank();

  for (int i=0; i<nprocs; ++i) {
    if ( comm_size != nprocs +1 )
      proc_reference_list_.h_view(i) = i;
    else {
      proc_reference_list_.h_view(i) = i;
      if ( i >= rank)
        proc_reference_list_.h_view(i) = i+1;
    }
  }
  Kokkos::deep_copy(proc_reference_list_.d_view,proc_reference_list_.h_view);

  policy = TeamPolicy(nparts, 1);
  Kokkos::parallel_for(policy, FillMigrateList(migrate_list, nparts, nprocs));

  Kokkos::DualView<LO*> count("count", nprocs+1);

  Kokkos::parallel_reduce(nparts,   CountParticleByProcessor(proc_reference_list_.d_view, migrate_list), &count.d_view(0));
  Kokkos::deep_copy(count.h_view, count.d_view);
  Kokkos::View<LO*> index("TestFillMigrateList::Index", nprocs);

  Kokkos::View<Particle*> migrate_particles("MigrateParticles::MigrateParticles", nparts);
  Kokkos::View<LO*,HostSpace> offset("TestFillMigrateList::offset", nprocs+1);

  offset(0) = 0;
  for (int n=1; n<nprocs+1; ++n)
    offset(n) = offset(n-1)+count.h_view(n-1);
  Kokkos::deep_copy(count.d_view, offset);


  PackParticlesForMigration packer(parts, migrate_list, migrate_particles, count.d_view);
  return packer;
}
TEUCHOS_UNIT_TEST( Int, TestFillMigrateList ) {

  int nprocs=13;
  int nparts=123;

  PackParticlesForMigration packer = create_packer(nprocs, nparts);

  Kokkos::DualView<LO*> count("count", nprocs+1);
  Kokkos::DualView<LO*> offset("offset", nprocs+1);

  count.d_view = packer.get_index();
  Kokkos::deep_copy(offset.h_view, count.d_view);

  Kokkos::parallel_for(nparts, packer);
  count.d_view = packer.get_index();

  Kokkos::deep_copy(count.h_view, count.d_view);

  for (int n=1; n<nprocs+1; ++n)
    TEST_EQUALITY_CONST(count.h_view(n-1), offset.h_view(n));
}


struct TestMigrationSetup {

  TestMigrationSetup(DataWarehouse &data) : data_(data), particles_(data_.particles()) {
    rank_ = particles_.comm_->getRank();
    int size = particles_.comm_->getSize();
    std::set<int> neighbors;
    for (int i=0; i<size; ++i)
      if ( i != rank_ )
        neighbors.insert(i);
    particles_.set_neighboring_procs(neighbors);
    data_.particles_reference() = particles_;

  }
  KOKKOS_INLINE_FUNCTION
   void operator () (const member_type & dev) const {
    GO i = dev.league_rank();
    int dest_rank = particles_.ielement(i);
    if (particles_.ielement(i) != rank_)
      particles_.add_particles_to_migrate_list(dev, &i, &dest_rank, 1);
  }
  DataWarehouse &data_;
  ParticleList particles_;

  int rank_;
};



struct CheckMyParts {
  CheckMyParts(ParticleList parts): particles_(parts), rank_(comm->getRank()) { }
  typedef int value_type;
  void operator() (std::size_t i, value_type &lsum) const {
    if ( particles_.ielement(i) != rank_)
      lsum = 1;
    lsum = 0;
  }

  ParticleList particles_;
  int rank_;

};

TEUCHOS_UNIT_TEST( Int, TestMigration ) {
  comm->barrier();
  int size = comm->getSize();
  if (size == 1 ) return;
  Teuchos::RCP<Map> map = Teuchos::RCP<Map >(new Map(size, 0, comm));

  DataWarehouse data;
  data.create_particle_list(1024);
  ParticleList parts = data.particles();
  parts.set_comm(comm);
  parts.set_element_map(map);
  data.particles_reference() = parts;

  int nparts = 13;
  TeamPolicy policy(nparts, 1);
  Kokkos::parallel_for(policy, StoreMultipleParticle(parts, nparts, size));

  TestMigrationSetup migrate_setup(data);
  Kokkos::parallel_for(policy, migrate_setup);

  MigrateParticles migrater(data.particles(), data.mesh());
  migrater.execute();

  int total_parts;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &nparts, &total_parts);

  TEST_EQUALITY_CONST(total_parts, size*nparts);

  int error=0;
  Kokkos::parallel_reduce(parts.num_particles(), CheckMyParts(parts), error);
  TEST_EQUALITY_CONST(error, 0);

}

