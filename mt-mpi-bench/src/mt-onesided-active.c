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

#define ITERATIONS   4
#define DATA_ITEMS   8

double
get_seconds()
{
    struct timeval now;

    gettimeofday(&now, NULL);

    const double seconds = (double) now.tv_sec;

    const double usecs = (double) now.tv_usec;

    return seconds + (usecs * 1.0e-6);
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

    char *hostname = (char *) malloc(sizeof(char) * 256);

    int hostnameLen = 0;

    MPI_Get_processor_name(hostname, &hostnameLen);

    MPI_Barrier(MPI_COMM_WORLD);

#ifdef MTBENCH_THREADED
    if (0 == rank) {
        printf("=========================================================\n");
        printf("MPI Multithreaded One-Sided Communication Benchmark\n");
        printf("Using %d threads\n", threads);
        printf("=========================================================\n");
        printf("\n");

        fflush(stdout);
    }
#else
    if (0 == rank) {
        printf("=========================================================\n");
        printf("MPI Single Threaded One-Sided Communication Benchmark\n");
        printf("=========================================================\n");
        printf("\n");

        fflush(stdout);
    }
#endif

    if (0 == rank) {
        printf("Number of MPI Ranks: %d\n", size);
        printf("\n");
        fflush(stdout);
    }

    printf("Rank: %7d on processor %s\n", rank, hostname);
    free(hostname);

    fflush(stdout);

    MPI_Barrier(MPI_COMM_WORLD);

    if (size % 2 != 0) {
        if (0 == rank) {
            fprintf(stderr, "Error: must use an even number of MPI ranks.\n");
        }

        MPI_Abort(MPI_COMM_WORLD, -2);
    }

    const int half_rank = size / 2;

    // Allocate 32M double values
    const int local_data_count = 1024 * 1024 * DATA_ITEMS;

    double *local_data = NULL;

    MPI_Alloc_mem(sizeof(double) * local_data_count, MPI_INFO_NULL,
                  &local_data);

    if (rank < half_rank) {
#ifdef MTBENCH_THREADED
#pragma omp parallel for
#endif
        for (int i = 0; i < local_data_count; i++) {
            local_data[i] = 1.0;
        }
    } else {
#ifdef MTBENCH_THREADED
#pragma omp parallel for
#endif
        for (int i = 0; i < local_data_count; i++) {
            local_data[i] = 0.0;
        }
    }

    if (0 == rank) {
        printf("\n");
        printf("Ranks  < %d will perform Put operations\n", half_rank);
        printf("Ranks >= %d will wait for fence\n\n", half_rank);
        printf("Creating MPI window ...\n");
    }

    MPI_Win local_window;

    MPI_Win_create(local_data, local_data_count * sizeof(double), sizeof(double),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &local_window);
    MPI_Win_fence(0, local_window);

    if (0 == rank) {
#ifdef MTBENCH_THREADED
        printf("Performing threaded one-sided benchmark ...\n");
#else
        printf("Performing one-sided benchmark ...\n");
#endif
    }

    if (rank < half_rank) {
        int i = local_data_count;

        MPI_Barrier(MPI_COMM_WORLD);
        const double start = get_seconds();

        for (int i = 0; i < ITERATIONS; i++) {
#ifdef MTBENCH_THREADED
#pragma omp parallel for
#endif
            for (int k = 0; k < local_data_count; k++) {
                MPI_Put(&local_data[k], 1, MPI_DOUBLE,
                        rank + half_rank, k, 1, MPI_DOUBLE, local_window);
            }

            MPI_Win_fence(0, local_window);
        }

        const double end = get_seconds();

        MPI_Barrier(MPI_COMM_WORLD);

        if (0 == rank) {
            const double time_taken = end - start;

            const double us_time_taken = time_taken * 1000.0 * 1000.0;

            const double put_performed =
                (double) (local_data_count * ITERATIONS);
            const double mb_exchanged =
                (local_data_count * sizeof(double) *
                 ITERATIONS) / (1024.0 * 1024.0);

            double globalPuts = (double) local_data_count * ITERATIONS;

            globalPuts *= (double) half_rank;

            printf
                ("Global:     8B                                                            %18.4f Put/s\n",
                 globalPuts / time_taken);
            printf("Rank 0:     8B %15d %15.4f us %18.7f MB/s %18.4f Put/s\n",
                   local_data_count, us_time_taken / put_performed,
                   mb_exchanged / time_taken,
                   (local_data_count * ITERATIONS) / time_taken);
        }
    } else {
        MPI_Barrier(MPI_COMM_WORLD);

        for (int i = 0; i < ITERATIONS; i++) {
            MPI_Win_fence(0, local_window);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Win_free(&local_window);
    MPI_Barrier(MPI_COMM_WORLD);

    if (half_rank == rank) {
        printf("Performing data check on recv side ...\n");

        // Perform checks on the recv side of the Puts
        for (int i = 0; i < local_data_count; i++) {
            if (local_data[i] != 1.0) {
                fprintf(stderr,
                        "Error: value at index %d should be 1.0 value is %f\n",
                        i, local_data[i]);
                exit(-1);
            }
        }

        printf("Data check PASSED SUCCESSFULLY.\n");
        fflush(stdout);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /////////////////////////////////////////////////////////////////////

    if (0 == rank) {
        printf("\n");
        printf("=========================================================\n");
    }

    MPI_Finalize();
    return 0;
}
