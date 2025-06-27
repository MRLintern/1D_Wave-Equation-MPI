#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include "wave.h"

// Simple test for the update function
void test_update() {
    int id, p;
    int master_num_nodes = 10;
    int num_dt = 10;
    double dt = 0.001;

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    // Determine local nodes
    int local_num_nodes = master_num_nodes / p;
    if (id == p - 1) {
        local_num_nodes += master_num_nodes % p;
    }

    // Allocate local array
    double *u_local = (double *)malloc(local_num_nodes * sizeof(double));
    assert(u_local != NULL);

    // Call update and check output
    double *result = update(id, p, master_num_nodes, local_num_nodes, num_dt, dt, u_local);
    assert(result != NULL);

    /*
    // Basic validation: should not contain NaNs
    for (int i = 0; i < local_num_nodes; i++) {
        assert(result[i] == result[i]);  // NaN fails this
    }

    //free(result);
    */
    printf("On Processor %d: test_update passed.\n", id);
}

// Simple test for the collect function
void test_collect() {
    int id, p;
    int master_num_nodes = 10;
    int num_dt = 10;
    double dt = 0.001;

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    int local_num_nodes = master_num_nodes / p;
    if (id == p - 1) {
        local_num_nodes += master_num_nodes % p;
    }

    double *u_local = (double *)malloc(local_num_nodes * sizeof(double));
    assert(u_local != NULL);

    // Fill dummy values for collection
    for (int i = 0; i < local_num_nodes; i++) {
        u_local[i] = id + 0.1 * i;
    }

    collect(id, p, master_num_nodes, local_num_nodes, num_dt, dt, u_local);
    printf("On Processor %d: test_collect passed.\n", id);

    free(u_local);
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    test_update();
    test_collect();

    MPI_Finalize();
    return 0;
}