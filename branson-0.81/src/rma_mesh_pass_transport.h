//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rma_mesh_pass_transport.h
 * \author Alex Long
 * \date   March 2 2016
 * \brief  Transport routine using one sided messaging and mesh-passing DD
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef transport_rma_mesh_pass_h_
#define transport_rma_mesh_pass_h_

#include <algorithm>
#include <mpi.h>
#include <numeric>
#include <queue>
#include <vector>

#include "RNG.h"
#include "constants.h"
#include "decompose_photons.h"
#include "info.h"
#include "mesh.h"
#include "mesh_pass_transport.h"
#include "mesh_rma_manager.h"
#include "message_counter.h"
#include "mpi_types.h"
#include "sampling_functions.h"
#include "timer.h"

//! Transport photons from a source object using the mesh-passing algorithm
// and one-sided messaging to fulfill requests for mesh data
std::vector<Photon>
rma_mesh_pass_transport(Source &source, Mesh &mesh, IMC_State &imc_state,
                        const IMC_Parameters &imc_parameters,
                        RMA_Manager &rma_manager, Tally_Manager &tally_manager,
                        Message_Counter &mctr, std::vector<double> &rank_abs_E,
                        std::vector<double> &rank_track_E,
                        const MPI_Types &mpi_types, const Info &mpi_info) {
  using Constants::event_type;
  using std::queue;
  using std::vector;
  // events
  using Constants::CENSUS;
  using Constants::EXIT;
  using Constants::KILL;
  using Constants::WAIT;

  uint32_t n_local = source.get_n_photon();
  uint32_t n_local_sourced = 0;

  uint32_t cell_id;

  double census_E = 0.0;
  double exit_E = 0.0;
  double dt = imc_state.get_next_dt();      //! For making current photons
  double next_dt = imc_state.get_next_dt(); //! For census photons

  RNG *rng = imc_state.get_rng();
  Photon phtn;

  // timing
  Timer t_transport;
  Timer t_rebalance_census;
  t_transport.start_timer("timestep transport");

  bool new_data = false;       //! New data flag is initially false
  std::vector<Cell> new_cells; // New cells from completed RMA requests

  // number of particles to run between MPI communication
  const uint32_t batch_size = imc_parameters.get_batch_size();

  event_type event;
  uint32_t wait_list_size;

  // for tallying off rank data
  std::unordered_map<uint32_t, double> off_rank_abs_E;

  //--------------------------------------------------------------------------//
  // main loop over photons
  //--------------------------------------------------------------------------//
  vector<Photon> census_list; //! Local end of timestep census list
  queue<Photon> wait_list;    //! Photons waiting for mesh data

  while (n_local_sourced < n_local) {

    uint32_t n = batch_size;

    while (n && n_local_sourced < n_local) {

      phtn = source.get_photon(rng, dt);
      n_local_sourced++;

      // get start cell, this only change with cell crossing event
      cell_id = phtn.get_cell();

      // if mesh available, transport and process, otherwise put on the
      // waiting list
      if (mesh.mesh_available(cell_id)) {
        event = transport_photon_mesh_pass(phtn, mesh, rng, next_dt, exit_E,
                                           census_E, rank_abs_E, rank_track_E,
                                           off_rank_abs_E);
        cell_id = phtn.get_cell();
      } else
        event = WAIT;

      if (event == CENSUS)
        census_list.push_back(phtn);
      else if (event == WAIT) {
        rma_manager.request_cell_rma(phtn.get_grip(), mctr);
        wait_list.push(phtn);
      }
      n--;
    } // end batch transport

    // process off rank tally data don't force send
    bool force_send = false;
    tally_manager.process_off_rank_tallies(mctr, off_rank_abs_E, force_send);

    // process mesh requests
    new_cells = rma_manager.process_rma_mesh_requests(mctr);
    new_data = !new_cells.empty();
    if (new_data)
      mesh.add_non_local_mesh_cells(new_cells, rma_manager.get_n_new_cells());
    // if data was received, try to transport photons on waiting list
    if (new_data) {
      wait_list_size = wait_list.size();
      for (uint32_t wp = 0; wp < wait_list_size; ++wp) {
        phtn = wait_list.front();
        wait_list.pop();
        cell_id = phtn.get_cell();
        if (mesh.mesh_available(cell_id)) {
          event = transport_photon_mesh_pass(phtn, mesh, rng, next_dt, exit_E,
                                             census_E, rank_abs_E, rank_track_E,
                                             off_rank_abs_E);
          cell_id = phtn.get_cell();
        } else
          event = WAIT;

        if (event == CENSUS)
          census_list.push_back(phtn);
        else if (event == WAIT) {
          rma_manager.request_cell_rma(phtn.get_grip(), mctr);
          wait_list.push(phtn);
        }
      } // end wp in wait_list
    }   // end if no data

  } // end while (n_local_source < n_local)

  //--------------------------------------------------------------------------//
  // Main transport loop finished, transport photons waiting for data
  //--------------------------------------------------------------------------//
  wait_list_size = wait_list.size();
  while (!wait_list.empty()) {

    // process off rank tally data don't force send
    bool force_send = false;
    tally_manager.process_off_rank_tallies(mctr, off_rank_abs_E, force_send);

    // process mesh requests
    new_cells = rma_manager.process_rma_mesh_requests(mctr);
    new_data = !new_cells.empty();
    if (new_data)
      mesh.add_non_local_mesh_cells(new_cells, rma_manager.get_n_new_cells());

    // if new data received or there are no active mesh requests, try to
    // transport waiting list (it could be that there are no active memory
    // requests because the request queue was full at the time)
    if (new_data || rma_manager.no_active_requests()) {
      wait_list_size = wait_list.size();
      for (uint32_t wp = 0; wp < wait_list_size; ++wp) {
        phtn = wait_list.front();
        wait_list.pop();
        cell_id = phtn.get_cell();
        if (mesh.mesh_available(cell_id)) {
          event = transport_photon_mesh_pass(phtn, mesh, rng, next_dt, exit_E,
                                             census_E, rank_abs_E, rank_track_E,
                                             off_rank_abs_E);
          cell_id = phtn.get_cell();
        } else
          event = WAIT;

        if (event == CENSUS)
          census_list.push_back(phtn);
        else if (event == WAIT) {
          rma_manager.request_cell_rma(phtn.get_grip(), mctr);
          wait_list.push(phtn);
        }
      } // end wp in wait_list
    }   // end if new_data
  }     // end while wait_list not empty

  // process off rank tally data and force send (force send requires completion
  // of all tallies before continuing
  bool force_send = true;
  tally_manager.process_off_rank_tallies(mctr, off_rank_abs_E, force_send);

  // record time of transport work for this rank
  t_transport.stop_timer("timestep transport");

  MPI_Barrier(MPI_COMM_WORLD);

  // append remote tally to the current tally
  tally_manager.add_remote_tally(rank_abs_E);

  // set the preffered census size to 10% of the user photon number and comb
  uint64_t max_census_photons = 0.1 * imc_parameters.get_n_user_photon();
  comb_photons(census_list, max_census_photons, rng);

  // all ranks have now finished transport set diagnostic quantities
  imc_state.set_exit_E(exit_E);
  imc_state.set_post_census_E(census_E);
  imc_state.set_network_message_counts(mctr);
  imc_state.set_rank_transport_runtime(
      t_transport.get_time("timestep transport"));

  // send the off-rank census back to ranks that own the mesh its on and receive
  // census particles that are on your mesh

  t_rebalance_census.start_timer("timestep rebalance_census");
  vector<Photon> rebalanced_census =
      rebalance_raw_census(census_list, mesh, mpi_types);
  t_rebalance_census.stop_timer("timestep rebalance_census");

  imc_state.set_rank_rebalance_time(
      t_rebalance_census.get_time("timestep rebalance_census"));

  // only do this if off rank census was separate from on rank census
  // census_list.insert(census_list.end(), rebalanced_census.begin(),
  //  rebalanced_census.end());

  // sort on census vectors by cell ID (global ID)
  sort(census_list.begin(), census_list.end());

  // set post census size after sorting and merging
  imc_state.set_census_size(census_list.size());

  return census_list;
}

#endif // def transport_rma_mesh_pass_h_
//----------------------------------------------------------------------------//
// end of transport_mesh_pass_rma.h
//----------------------------------------------------------------------------//
