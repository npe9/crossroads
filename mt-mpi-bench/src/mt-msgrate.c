/*
<<<<<<< HEAD
 * Copyright 2006 - 2016 Sandia Corporation.
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia
 * Corporation, the U.S. Government retains certain rights in
 * this software.
=======
 * Copyright 2015-2016 Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
 * Government retains certain rights in this software.
>>>>>>> 92613e83fe31d76f883007369f2738c71a96fbc0
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
<<<<<<< HEAD
 *
=======
>>>>>>> 92613e83fe31d76f883007369f2738c71a96fbc0
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include <mpi.h>

#ifdef MTBENCH_THREADED
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

#define ITERATIONS   5000
#define MAX_MSG_SIZE 16384

double
get_seconds()
{
    struct timeval now;

    gettimeofday(&now, NULL);

    const double seconds = (double) now.tv_sec;
    const double usecs = (double) now.tv_usec;

    return seconds + (usecs * 1.0e-6);
}

void
perform_recv(const int source, const int world_size,
             double *buffer, const int msg_size, MPI_Request * req)
{

    int recv_err = MPI_Irecv(buffer, msg_size, MPI_DOUBLE, source, 0,
                             MPI_COMM_WORLD, req);

    if (MPI_SUCCESS != recv_err) {
        MPI_Abort(MPI_COMM_WORLD, -4);
    }
}

void
perform_send(const int dest, const int world_size,
             double *buffer, const int msg_size, MPI_Request * req)
{

    int send_err =
        MPI_Isend(buffer, msg_size, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD,
                  req);

    if (MPI_SUCCESS != send_err) {
        MPI_Abort(MPI_COMM_WORLD, -4);
    }
}

int
main(int argc, char *argv[])
{

#ifdef MTBENCH_THREADED
    int mpi_provided = MPI_THREAD_SINGLE;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpi_provided);

    if (MPI_THREAD_MULTIPLE != mpi_provided) {
        fprintf(stderr, "Error: Requested MPI_THREAD_MULTIPLE but was given ");

        switch (mpi_provided) {
        case MPI_THREAD_SINGLE:
            fprintf(stderr, "MPI_THREAD_SINGLE\n");
            break;
        case MPI_THREAD_FUNNELED:
            fprintf(stderr, "MPI_THREAD_FUNNELED\n");
            break;
        case MPI_THREAD_SERIALIZED:
            fprintf(stderr, "MPI_THREAD_SERIALIZED\n");
            break;
        }

        return -1;
    }
#else
    MPI_Init(&argc, &argv);
#endif

    /////////////////////////////////////////////////////////////////////

#if defined(MTBENCH_THREADED) && defined(_OPENMP)
    const int threads = omp_get_max_threads();
#else
    const int threads = 1;
#endif

    int rank = 0;
    int size = 0;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#ifdef MTBENCH_THREADED
    if (0 == rank) {
        printf("=========================================================\n");
        printf("MPI Multithreaded Communication Benchmark\n");
        printf("Using %d threads\n", threads);
        printf("=========================================================\n");
        printf("\n");
    }
#else
    if (0 == rank) {
        printf("=========================================================\n");
        printf("MPI Single Threaded Communication Benchmark\n");
        printf("=========================================================\n");
        printf("\n");
    }
#endif

    const int half_rank = size / 2;

    if (0 == rank) {
        printf("Number of MPI Ranks: %d\n", size);
        printf("Each rank requires %6.3f MB for MPI buffers.\n",
               (MAX_MSG_SIZE * ITERATIONS * sizeof(double) /
                (1024.0 * 1024.0)));
        printf("\n");
        printf("Per-Rank Performance Data (measured from rank 0):\n\n");
        printf("%10s %15s    %18s      %15s        %10s\n",
               "Size", "Latency", "Local Bandwidth", "Local Msg Rate",
               "Glbl. MRate");
    }

    if (size % 2 > 0) {
        if (0 == rank) {
            fprintf(stderr, "Error: must use an even number of MPI ranks.\n");
        }

        MPI_Abort(MPI_COMM_WORLD, -2);
    }
#ifdef MTBENCH_THREADED
#pragma omp parallel
    {
#endif
        if (rank < half_rank) {
            for (int i = 1; i <= MAX_MSG_SIZE; i = i * 2) {
                double *buffer =
                    (double *) malloc(sizeof(double) * i * ITERATIONS);

                for (int k = 0; k < ITERATIONS; k++) {
                    for (int j = 0; j < i; j++) {
                        buffer[(k * i) + j] = (double) k *i;
                    }
                }

                MPI_Request *requests = (MPI_Request *)
                    malloc(sizeof(MPI_Request) * ITERATIONS);

                MPI_Status *statuses = (MPI_Status *) malloc(sizeof(MPI_Status)
                                                             * ITERATIONS);

#ifdef MTBENCH_THREADED
#pragma omp barrier
#pragma omp single
                {
#endif
                    MPI_Barrier(MPI_COMM_WORLD);
#ifdef MTBENCH_THREADED
                }
#endif
                const double start = get_seconds();

                for (int k = 0; k < ITERATIONS; k++) {
                    perform_recv(rank + half_rank, size,
                                 &buffer[k * i], i, &requests[k]);
                }

                MPI_Waitall(ITERATIONS, requests, statuses);

#ifdef MTBENCH_THREADED
#pragma omp barrier
#pragma omp single
                {
#endif
                    MPI_Barrier(MPI_COMM_WORLD);
#ifdef MTBENCH_THREADED
                }
#endif
                const double end = get_seconds();

                if (0 == rank) {
#ifdef MTBENCH_THREADED
#pragma omp master
                    {
#endif
                        const double time_taken = end - start;

                        const double us_time_taken_itr =
                            (time_taken /
                             (double) ITERATIONS) * 1000.0 * 1000.0;
                        const double mb_exchanged = ((double)
                                                     (i * ITERATIONS) *
                                                     2.0 * sizeof(double) *
                                                     threads) / (1024.0 *
                                                                 1024.0);
                        const double msgs_xchng = ITERATIONS * threads;

                        const double msgs_xchng_global =
                            (((size / 2) * msgs_xchng) / 1000000.0);

                        printf
                            ("%10d %15.4f us %18.7f MB/s %15.4f msg/s %12.4f Mmsgs/s\n",
                             i, us_time_taken_itr, mb_exchanged / time_taken,
                             msgs_xchng / time_taken,
                             msgs_xchng_global / time_taken);
#ifdef MTBENCH_THREADED
                    }
#endif
                }

                free(buffer);
                free(requests);
                free(statuses);
            }
        } else {
            for (int i = 1; i <= MAX_MSG_SIZE; i = i * 2) {

                MPI_Request *requests = (MPI_Request *)
                    malloc(sizeof(MPI_Request) * ITERATIONS);

                MPI_Status *statuses = (MPI_Status *) malloc(sizeof(MPI_Status)
                                                             * ITERATIONS);

                double *buffer =
                    (double *) malloc(sizeof(double) * i * ITERATIONS);

                for (int k = 0; k < ITERATIONS; k++) {
                    for (int j = 0; j < i; j++) {
                        buffer[(k * i) + j] = 0;
                    }
                }

#ifdef MTBENCH_THREADED
#pragma omp barrier
#pragma omp single
                {
#endif
                    MPI_Barrier(MPI_COMM_WORLD);
#ifdef MTBENCH_THREADED
                }
#endif
                for (int k = 0; k < ITERATIONS; k++) {
                    perform_send(rank - half_rank, size,
                                 &buffer[k * i], i, &requests[k]);
                }

                MPI_Waitall(ITERATIONS, requests, statuses);

#ifdef MTBENCH_THREADED
#pragma omp barrier
#pragma omp single
                {
#endif
                    MPI_Barrier(MPI_COMM_WORLD);
#ifdef MTBENCH_THREADED
                }
#endif

                free(buffer);
                free(requests);
                free(statuses);
            }
        }
#ifdef MTBENCH_THREADED
    }
#endif

    /////////////////////////////////////////////////////////////////////

    if (0 == rank) {
        printf("\n");
        printf("=========================================================\n");
    }

    MPI_Finalize();
    return 0;
}
