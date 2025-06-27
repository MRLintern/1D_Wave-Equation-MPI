#ifndef WAVE_H
#define WAVE_H

// Function to update wave field for each MPI process
double *update(int id, int p, int master_num_nodes, int local_num_nodes, int num_dt, double dt, double *u_local);

// Function to collect final results on master process
void collect(int id, int p, int master_num_nodes, int local_num_nodes, int num_dt, double dt, double *u_local);

// Exact initial condition function (user-defined elsewhere)
double exact(double x, double t);

// Function to print a timestamp
void timestamp(void);

#endif