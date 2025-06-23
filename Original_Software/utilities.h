#ifndef UTILITIES
#define UTILITIES

//--function declarations--//

//function to advance/move solution to the next time-step
double* update(int id, int p, int master_num_nodes, int local_num_nodes, int num_dt, double dt);

//function to send calculated data from the local process to the master process
void collect(int id, int p, int master_num_nodes, int local_num_nodes, int num_dt, double dt, double u_local[]);

//for evaluating time derivative of the solution
double dudt(double x, double t);

//function to evaluate exact solution
double exact(double x, double t);

//function for time considerations
void timestamp();

#endif
