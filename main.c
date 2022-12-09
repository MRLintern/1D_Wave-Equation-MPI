/*

Program solves the 2D Wave Equation.
Equation is discretized using a Finite Difference Stencil.
The resulting System of Algebraic Equations is solved
using the Jacobi Iterative Method which is parallelized using MPI.

The program was run with a different number of processors
along with different compiler optimization flags: -O, -O0, -O1, -O2 & -O3.

The application is run by using the shell script, wave.sh
The compiler opt. flag and the number of processors (-np)
can be changed in the shell script.

Below is the finite difference equation, initial and boundary conditions.

I.C.:

u(x, 0) = I(x); shape of the wave
d/dt (x, 0) = 0; initial velocity is 0.

B.C.:

Dirichlet B.C.; the wave is fixed at each end

u(L, t) = 0; L is the length of the spacial domain
u(0, t) = 0; 

Equation:

The stability condition for the solution: Courant–Friedrichs–Lewy condition (CFL)

C = c*dt/dx. Note: in the code, C = CFL.

We assume that we know the solution at the starting node. However, we don't know the
solution at next node. This is what we're trying to find.

Note: n represents the nodes in the time direction. i represents the nodes in the space direction.

u(i,n+1) = -u(i,n-1) + 2u(i,n) + C^2(u(i+1,n) - 2u(i,n) + u(i-1,n))

*/

//libraries
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<mpi.h>

//--function declarations--//
////////////////////////////

double* update(int id, int p, int n_global, int n_local, int nsteps, double dt);
void collect(int id, int p, int n_global, int n_local, int nsteps, double dt, double u_local[]);
double dudt(double x, double t);
double exact(double x, double t);
void timestamp();

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/

{
  double dt = 0.00125;
  int i_global_hi;
  int i_global_lo;
  int id;
  int n_global = 401;
  int n_local;
  int nsteps = 4000;
  int p;
  double *u1_local;
  double wtime;
/* 
  Initialize MPI.
*/
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  MPI_Comm_size(MPI_COMM_WORLD, &p);

  if (id == 0) 
  {
    timestamp();
    printf("Using %d processes.\n", p);
    printf("Using a total of %d points.\n", n_global);
    printf("Using %d time steps of size %g.\n", nsteps, dt);
  }

  wtime = MPI_Wtime();
/*
  Determine N_LOCAL.
*/
  i_global_lo = (id*(n_global - 1))/p;
  i_global_hi = ((id + 1)*(n_global - 1))/p;
  
  if (0 < id)
  {
    i_global_lo = i_global_lo - 1;
  }
  n_local = i_global_hi + 1 - i_global_lo;
/* 
  Update N_LOCAL values.
*/
  u1_local = update(id, p, n_global, n_local, nsteps, dt);
/* 
  Collect local values into global array.
*/
  collect(id, p, n_global, n_local, nsteps, dt, u1_local);
/*
  Report elapsed wallclock time.
*/
  wtime = MPI_Wtime() - wtime;
  if (id == 0)
  {
    printf("\n");
    printf("Elapsed wallclock time was %g seconds\n", wtime);
  }
/*
  Terminate MPI.
*/
  MPI_Finalize();
/*
  Free memory.
*/
  free(u1_local);
/*
  Terminate.
*/
  if (id == 0)
  {
    printf("\n");
    printf("Normal end of execution.\n");
    timestamp();
  }

  return 0;
}
/******************************************************************************/

double *update(int id, int p, int n_global, int n_local, int nsteps, double dt) 

/******************************************************************************/
/*
  Purpose:

    UPDATE advances the solution a given number of time steps.

*/
{
  double alpha;
  double c;
  double dx;
  int i;
  int i_global;
  int i_global_hi;
  int i_global_lo;
  int i_local;
  int i_local_hi;
  int i_local_lo;
  int ltor = 20;
  int rtol = 10;
  MPI_Status status;
  double t;
  double *u0_local;
  double *u1_local;
  double *u2_local;
  double x;
/*
  Determine the value of ALPHA.
*/
  c = 1.0;
  dx = 1.0/(double)(n_global - 1);
  alpha = c*dt/dx;

  if (1.0 <= fabs(alpha))
  {
    if (id == 0)
    {
      fprintf(stderr, "\n");
      fprintf(stderr, "UPDATE - Warning!\n");
      fprintf(stderr, "  1 <= |ALPHA| = | C * dT / dX |.\n");
      fprintf(stderr, "  C = %g\n", c);
      fprintf(stderr, "  dT = %g\n", dt);
      fprintf(stderr, "  dX = %g\n", dx);
      fprintf(stderr, "  ALPHA = %g\n", alpha);
      fprintf(stderr, "  Computation will not be stable!\n");
    }
    MPI_Finalize();
    exit(1);
  }
/*
  The global array of N_GLOBAL points must be divided up among the processes.
  Each process stores about 1/P of the total + 2 extra slots.
*/
  i_global_lo = (id*(n_global - 1))/p;
  i_global_hi = (( id + 1)*(n_global - 1))/p;
  
  if (0 < id)
  {
    i_global_lo = i_global_lo - 1;
  }

  i_local_lo = 0;
  i_local_hi = i_global_hi - i_global_lo;

  u0_local = (double*)malloc(n_local*sizeof(double));
  u1_local = (double*)malloc(n_local*sizeof(double));
  u2_local = (double*)malloc(n_local*sizeof(double));

  t = 0.0;
  
  for (i_global = i_global_lo; i_global <= i_global_hi; i_global++) 
  {
    x = (double)(i_global)/(double)(n_global - 1);
    i_local = i_global - i_global_lo;
    u1_local[i_local] = exact(x, t);
  }

  for (i_local = i_local_lo; i_local <= i_local_hi; i_local++)
  {
    u0_local[i_local] = u1_local[i_local];
  }
/* 
  Take NSTEPS time steps.
*/
  for (i = 1; i <= nsteps; i++)
  {
    t = dt*(double)i;
/* 
  For the first time step, we need to use the initial derivative information.
*/
    if (i == 1)
    {
      for (i_local = i_local_lo + 1; i_local < i_local_hi; i_local++) 
      {
        i_global = i_global_lo + i_local;
        x = (double)(i_global)/(double)(n_global - 1);
        u2_local[i_local] = 
          
          + 0.5*alpha*alpha*u1_local[i_local-1]
          + (1.0 - alpha*alpha)*u1_local[i_local] 
          +  0.5*alpha*alpha*u1_local[i_local+1]
          +  dt*dudt(x, t);
      }
    }
/* 
  After the first time step, we can use the previous two solution estimates.
*/
    else
    {
      for (i_local = i_local_lo + 1; i_local < i_local_hi; i_local++) 
      {
        u2_local[i_local] =
          
          + alpha*alpha*u1_local[i_local-1]
          + 2.0*(1.0 - alpha*alpha)*u1_local[i_local] 
          + alpha*alpha*u1_local[i_local+1]
          - u0_local[i_local];
      }
    }
/* 
  Exchange data with "left-hand" neighbor. 
*/
    if (0 < id) 
    {
      MPI_Send(&u2_local[i_local_lo+1], 1, MPI_DOUBLE, id - 1, rtol, MPI_COMM_WORLD);
      MPI_Recv(&u2_local[i_local_lo], 1, MPI_DOUBLE, id - 1, ltor, MPI_COMM_WORLD, &status );
    }
    else
    {
      x = 0.0;
      u2_local[i_local_lo] = exact(x, t);
    }
/* 
  Exchange data with "right-hand" neighbor.
*/
    if (id < p - 1) 
    {
      MPI_Send(&u2_local[i_local_hi-1], 1, MPI_DOUBLE, id + 1, ltor, MPI_COMM_WORLD );
      MPI_Recv(&u2_local[i_local_hi], 1, MPI_DOUBLE, id + 1, rtol, MPI_COMM_WORLD, &status);
    }
    else
    {
      x = 1.0;
      u2_local[i_local_hi] = exact(x, t);
    }
/*
  Shift data for next time step.
*/
    for (i_local = i_local_lo; i_local <= i_local_hi; i_local++)
    {
      u0_local[i_local] = u1_local[i_local];
      u1_local[i_local] = u2_local[i_local];
    }
  }
/*
  Free memory.
*/
  free(u0_local);
  free(u2_local);

  return u1_local;
}
/******************************************************************************/

void collect ( int id, int p, int n_global, int n_local, int nsteps, double dt, double u_local[]) 

/******************************************************************************/
/*
  Purpose:

    COLLECT has workers send results to the master, which prints them.

  
*/
{
  int buffer[2];
  int collect1 = 10;
  int collect2 = 20;
  int i;
  int i_global;
  int i_global_hi;
  int i_global_lo;
  int i_local;
  int i_local_hi;
  int i_local_lo;
  int n_local2;
  MPI_Status status;
  double t;
  double *u_global;
  double x;

  i_global_lo = (id*(n_global - 1))/p;
  i_global_hi = ((id + 1)*(n_global - 1)/p;
                 
  if (0 < id)
  {
    i_global_lo = i_global_lo - 1;
  }

  i_local_lo = 0;
  i_local_hi = i_global_hi - i_global_lo;
/* 
  Master collects worker results into the U_GLOBAL array.
*/
  if ( id == 0 )
  {
/*
  Create the global array.
*/
    u_global = (double*)malloc(n_global*sizeof(double));
/*
  Copy the master's results into the global array.
*/
    for (i_local = i_local_lo; i_local <= i_local_hi; i_local++)
    {
      i_global = i_global_lo + i_local - i_local_lo;
      u_global[i_global] = u_local[i_local];
    }
/*
  Contact each worker.
*/
    for (i = 1; i < p; i++) 
    {
/*
  Message "collect1" contains the global index and number of values.
*/
      MPI_Recv(buffer, 2, MPI_INT, i, collect1, MPI_COMM_WORLD, &status);
      i_global_lo = buffer[0];
      n_local2 = buffer[1];

      if (i_global_lo < 0)
      {
        fprintf(stderr, "Illegal I_GLOBAL_LO = %d\n", i_global_lo);
        exit ( 1 );
      }
      else if (n_global <= i_global_lo + n_local2 - 1)
      {
        fprintf(stderr, "  Illegal I_GLOBAL_LO + N_LOCAL2 = %d\n", i_global_lo + n_local2);
        exit ( 1 );
      }
/*
  Message "collect2" contains the values.
*/
      MPI_Recv(&u_global[i_global_lo], n_local2, MPI_DOUBLE, i, collect2, MPI_COMM_WORLD, &status);
    }
/*
  Print the results.
*/
    t = dt*(double)nsteps;
    printf("\n");
    printf("    I      X     F(X)   Exact\n");
    printf("\n");
    for (i_global = 0; i_global < n_global; i_global++) 
    {
      x = (double)(i_global)/(double)(n_global - 1);
      printf("  %3d  %6.3f  %6.3f  %6.3f\n", i_global, x, u_global[i_global], exact (x, t));
    }

    free(u_global);
  }
/*
  Workers send results to process 0.
*/
  else
  {
/*
  Message "collect1" contains the global index and number of values.
*/
    buffer[0] = i_global_lo;
    buffer[1] = n_local;
    MPI_Send(buffer, 2, MPI_INT, 0, collect1, MPI_COMM_WORLD);
/*
  Message "collect2" contains the values.
*/
    MPI_Send(u_local, n_local, MPI_DOUBLE, 0, collect2, MPI_COMM_WORLD);
  }

  return;
}
/******************************************************************************/

double exact ( double x, double t )

/******************************************************************************/
/*
  Purpose:

    EXACT evaluates the exact solution

 
*/
{
  const double c = 1.0;
  const double pi = 3.141592653589793;
  double value;

  value = sin(2.0*pi*(x - c*t));

  return value;
}
/******************************************************************************/

double dudt(double x, double t)

/******************************************************************************/
/*
  Purpose:

    DUDT evaluates the partial derivative dudt.

  
*/
{
  const double c = 1.0;
  const double pi = 3.141592653589793;
  double value;

  value = - 2.0*pi*c*cos(2.0*pi*(x - c*t));

  return value;
}
/******************************************************************************/

void timestamp()

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

 
*/
{
#define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm* tm;
  time_t now;

  now = time(NULL);
  tm = localtime(&now);

  strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

  printf ("%s\n", time_buffer);

  return;
  
#undef TIME_SIZE
}

