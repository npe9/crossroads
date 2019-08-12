!-----------------------------------------------------------------------
!
! MODULE: corsweep_module
!> @brief
!> This module controls the mesh sweep scheduling. It directs the flow
!> of KBA stages for different octants, groups, spatial work chunks.
!
!-----------------------------------------------------------------------

MODULE corsweep_module

  USE global_module, ONLY: i_knd, zero, one, l_knd

  USE data_module, ONLY: ng

  USE geom_module, ONLY: nc, ny, nz, jdim, kdim

  USE sn_module, ONLY: nang

  USE control_module, ONLY: inrdone, ncor, nops, yzstg,                &
    corner_loop_order

  USE octsweep_module, ONLY: octsweep

  USE solvar_module, ONLY: flkx, flky, flkz, fmin, fmax, flux0, fluxm, &
    jb_in, kb_in, jb_out, kb_out

  USE plib_module, ONLY: nthreads, waitinit, ichunk, firsty, lasty,    &
    firstz, lastz, iproc, npey, npez

  USE thrd_comm_module, ONLY: corsweep_recv_bdry, corsweep_test_pick,  &
    corsweep_send_bdry, assign_thrd_set

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: corsweep

  SAVE


  CONTAINS


  SUBROUTINE corsweep ( t, do_grp, ng_per_thrd, nnstd_used, grp_act )

!-----------------------------------------------------------------------
!
! Driver for the mesh sweeps. Manages the loops over octant pairs.
!
!-----------------------------------------------------------------------

    INTEGER(i_knd), INTENT(IN) :: t

    INTEGER(i_knd), INTENT(INOUT) :: ng_per_thrd, nnstd_used

    INTEGER(i_knd), DIMENSION(ng), INTENT(INOUT) :: do_grp

    INTEGER(i_knd), DIMENSION(ng,nthreads), INTENT(INOUT) :: grp_act
!_______________________________________________________________________
!
!   Local variables
!_______________________________________________________________________

    INTEGER(i_knd) :: jd, kd, g, n, nulreq, prev_c(2), cor, i,         &
      dummy(2)=(/0, 0/), nc_tot, nstages

    INTEGER(i_knd), DIMENSION(ncor) :: mtag, ic_op, gc, g_op

    INTEGER(i_knd), DIMENSION(2*ncor) :: reqr, reqs

    LOGICAL(l_knd), DIMENSION(ncor) :: done, callrecv
!_______________________________________________________________________
!
!   Assign the work to threads. Main level threads always applied to
!   energy groups. Apply nested threads additionally to groups if
!   swp_typ is 0. Apply nested threads to mini-KBA if swp_typ is 1.
!_______________________________________________________________________

  !$OMP MASTER

    do_grp = 1
    WHERE ( inrdone ) do_grp = 0

    CALL assign_thrd_set ( do_grp, ng, ng_per_thrd, 0, nnstd_used,     &
      grp_act )

  !$OMP END MASTER
  !$OMP BARRIER
!_______________________________________________________________________
!
!   Initialize the request arrays for receive (reqr) and send (reqs)
!_______________________________________________________________________

    CALL waitinit ( reqr, SIZE( reqr ) )
    CALL waitinit ( reqs, SIZE( reqs ) )
    nulreq = reqs(1)
!_______________________________________________________________________
!
!   Clean up and initialize some values.
!_______________________________________________________________________

    nc_tot = 0

    clean_loop: DO n = 1, ng_per_thrd

      g = grp_act(n,t)
      IF ( g == 0 ) EXIT clean_loop

      fmin(g) = HUGE( one )
      fmax(g) = zero

      flkx(:,:,:,g) = zero
      flky(:,:,:,g) = zero
      flkz(:,:,:,g) = zero

      flux0(:,:,:,g) = zero
      fluxm(:,:,:,:,g) = zero

      nc_tot = nc_tot + nc

    END DO clean_loop

    nc_tot = 2*nc_tot
    nstages = nc_tot + npey + npez - 2
!_______________________________________________________________________
!
!   Initialize all the counters that will control the concurrent octant
!   sweeps.
!
!   done(ncor) - flag for when all groups done for a starting corner
!   callrecv(ncor) - which starting corner needs to call for another msg
!   mtag(ncor) - message tag for each starting corner
!   prev_c(2) - what was the previous corner this process did
!   ic_op(ncor) - chunk counter for each corner over octant-pair
!   gc(ncor) - group corresponding to each starting corner
!   g_op(ncor) - counter for each corner over octant-pair and groups
!
!_______________________________________________________________________

    done     = .FALSE.
    callrecv = .TRUE.

    mtag   = 0
    prev_c = 0

    ic_op = 1
    gc    = 1
    g_op  = 1
!_______________________________________________________________________
!
!   Are there more threads than remaining groups to be swept.
!_______________________________________________________________________

    IF ( grp_act(1,t) == 0 ) done = .TRUE.
!_______________________________________________________________________
!
!   Loop over all chunks: all groups, all octants. Simultaneous sweeps
!   from starting corners (2 in 2D, 4 in 3D). Exit when all work done.
!_______________________________________________________________________

    done_loop: DO

      IF ( ALL( done ) ) EXIT
!_______________________________________________________________________
!
!     Call to non-blocking receive for the corner that has
!     callrecv=true: first pass that is all starting corners; subsequent
!     passes only one corner is true.
!_______________________________________________________________________

      CALL corsweep_recv_bdry ( jdim, kdim, ncor, ng, ic_op, gc,       &
        grp_act(:,t), callrecv, mtag, reqr, SIZE( reqr ), nang, ichunk,&
        ny, nz, jb_in, kb_in )
!_______________________________________________________________________
!
!     Determine next corner to be swept according to available message,
!     safe buffer space, and remaining stages for given direction.
!_______________________________________________________________________

      CALL corsweep_test_pick ( ncor, corner_loop_order, yzstg,        &
        nstages, nulreq, prev_c, g_op, gc, ng, grp_act(:,t), reqr,     &
        SIZE( reqr ), done, reqs, SIZE( reqs ), cor )
!_______________________________________________________________________
!
!     Set up the sweep for this corner.
!_______________________________________________________________________

      g = grp_act(gc(cor),t)

      SELECT CASE ( cor )
        CASE ( 1 )
          jd = 1; kd = 1
        CASE ( 2 )
          jd = 2; kd = 1
        CASE ( 3 )
          jd = 1; kd = 2
        CASE ( 4 )
          jd = 2; kd = 2
      END SELECT
!_______________________________________________________________________
!
!     Sweep the task defined by the group, chunk, octant. Send the
!     outgoing data downstream.
!_______________________________________________________________________

      CALL octsweep ( g, jd, kd, 0, dummy, 2, ic_op(cor), cor )

      CALL corsweep_send_bdry ( jd, kd, nang, ichunk, ny, nz, ncor,    &
        mtag(cor), reqs(cor), reqs(ncor+cor),                          &
        jb_out(:,:,:,cor,ic_op(cor),g), kb_out(:,:,:,cor,ic_op(cor),g) )
!_______________________________________________________________________
!
!     Increment the counters, and reset control variables.
!_______________________________________________________________________

!      IF ( prev_c(1) == cor ) THEN
!        prev_c(2) = prev_c(2) + 1
!      ELSE
!        prev_c(2) = 1
!      END IF
!      prev_c(1) = cor

      mtag(cor) = 0

      IF ( ic_op(cor) /= 2*nc ) THEN
        ic_op(cor) = ic_op(cor) + 1
      ELSE
        ic_op(cor) = 1
        gc(cor) = gc(cor) + 1
      END IF

      IF ( g_op(cor) /= nc_tot ) THEN
        g_op(cor) = g_op(cor) + 1
        callrecv(cor) = .TRUE.
      ELSE
        done(cor) = .TRUE.
      END IF
!_______________________________________________________________________
!
!   End loop over all operations
!_______________________________________________________________________

    END DO done_loop
!_______________________________________________________________________
!_______________________________________________________________________

  END SUBROUTINE corsweep


END MODULE corsweep_module
