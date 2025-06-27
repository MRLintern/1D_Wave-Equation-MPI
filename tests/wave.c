#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <time.h> 
#include "wave.h"


// Print a timestamp to indicate when the program started or finished.
void timestamp(void) {

    time_t now;
    time(&now);

    printf("%s", ctime(&now));  // Print formatted timestamp
}

// exact solution to wave equation: f(x,t) = sin(2.0*pi*(x - t))
double exact(double x, double t) {

    return sin(2.0 * M_PI * (x - t));
}

/*

 * Update the wave field using a finite difference method.
 * This function simulates the 1D wave equation u_tt = c^2 * u_xx
 * using a simple second-order finite difference in time and space.
 *
 * param id                 MPI rank of the current process.
 * param p                  Total number of MPI processes.
 * param master_num_nodes  Total number of spatial nodes in the full domain.
 * param local_num_nodes   Number of nodes this process is responsible for.
 * param num_dt            Number of time steps to simulate.
 * param dt                Time step size.
 * param u_local           Pointer to the local solution array.
 * return                  Pointer to updated local solution array.
 
 */

double *update(int id, int p, int master_num_nodes, int local_num_nodes, int num_dt, double dt, double *u_local) {

    double dx = 1.0 / (master_num_nodes - 1);   // Spatial step size
    double wave_speed = 1.0;                    // Speed of wave propagation
    double CFL = dt / dx;                       // Courant–Friedrichs–Lewy number

    // Ensure the CFL condition for stability is met
    if (CFL * wave_speed >= 1.0) {

        if (id == 0) {

            fprintf(stderr, "Error: CFL condition violated. CFL = %f\n", CFL * wave_speed);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Allocate arrays for the previous (u_old), current (u), and next (u_new) time steps
    double *u_old = (double *)calloc(local_num_nodes, sizeof(double));
    double *u = (double *)calloc(local_num_nodes, sizeof(double));
    double *u_new = (double *)calloc(local_num_nodes, sizeof(double));

    // Set initial condition: u(x, 0) = exact(x, 0)
    for (int i = 0; i < local_num_nodes; i++) {

        int global_index = id * (master_num_nodes - 1) / p + i;

        double x = global_index * dx;

        u[i] = exact(x, 0.0);

        u_old[i] = u[i];  // Assume zero initial velocity: u_t(x, 0) = 0
    }

    // Time-stepping loop
    for (int t = 0; t < num_dt; t++) {

        // Exchange ghost cells with neighboring processes
        if (id < p - 1) {

            MPI_Send(&u[local_num_nodes - 2], 1, MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&u[local_num_nodes - 1], 1, MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if (id > 0) {

            MPI_Recv(&u[0], 1, MPI_DOUBLE, id - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&u[1], 1, MPI_DOUBLE, id - 1, 0, MPI_COMM_WORLD);
        }

        // Update the interior points using finite difference scheme
        for (int i = 1; i < local_num_nodes - 1; i++) {

            u_new[i] = 2 * u[i] - u_old[i] + (CFL * wave_speed) * (CFL * wave_speed) * (u[i - 1] - 2 * u[i] + u[i + 1]);
        }

        // Apply boundary conditions (Dirichlet: u=0)
        if (id == 0) { u_new[1] = 0.0; }
        if (id == p - 1) { u_new[local_num_nodes - 2] = 0.0; }

        // Rotate time levels
        double *temp = u_old;

        u_old = u;
        u = u_new;
        u_new = temp;
    }

    // Copy final result into u_local
    memcpy(u_local, u, local_num_nodes * sizeof(double));

    // Free intermediate buffers
    free(u_old);
    free(u_new);

    return u_local;
}

/*

 * Gather the final results from all processes
 *
 * param id                 MPI rank of current process.
 * param p                  Total number of processes.
 * param master_num_nodes  Total number of nodes across the domain.
 * param local_num_nodes   Number of nodes handled by each process.
 * param num_dt            Total number of time steps.
 * param dt                Time step size.
 * param u_local           Pointer to this process’s final solution data.

 */
void collect(int id, int p, int master_num_nodes, int local_num_nodes, int num_dt, double dt, double *u_local) {

    // Rank 0 process will collect and save the entire domain solution
    if (id == 0) {

        double *u_global = (double *)malloc(master_num_nodes * sizeof(double));

        double dx = 1.0 / (master_num_nodes - 1);

        // Copy local data (excluding ghost cells)
        for (int i = 0; i < local_num_nodes - 1; i++) {

            u_global[i] = u_local[i];
        }

        int offset = local_num_nodes - 1;

        // Receive from other processes
        for (int source = 1; source < p; source++) {

            int recv_count = ((source + 1) * (master_num_nodes - 1)) / p - (source * (master_num_nodes - 1)) / p;

            double *recv_buf = (double *)malloc(recv_count * sizeof(double));

            MPI_Recv(recv_buf, recv_count, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for (int i = 0; i < recv_count; i++) {
                u_global[offset + i] = recv_buf[i];
            }

            offset += recv_count;
            free(recv_buf);
        }

        // Create x array for CSV output; NO CSV FILE CREATED OR WRITTEN TO; THIS IS A TEST CASE
        double *x = (double *)malloc(master_num_nodes * sizeof(double));

        for (int i = 0; i < master_num_nodes; i++) {

            x[i] = i * dx;
        }

        free(x);
        free(u_global);

    } else {

        // Non-master processes send data back to rank 0
        int send_count = local_num_nodes - 1;

        MPI_Send(u_local, send_count, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
}
