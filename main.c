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

//enter the main function
int main (int argc, char *argv[])
{
  
  //time step
  double dt = 0.00125;
  
  //--master process arrays
  int i_master_max; //max value
  int i_master_min; //min value
  
  //MPI process I.D.
  int id;
  
  //number of nodes for master process
  int master_num_nodes = 401;
  
  //number of nodes for local process
  int local_num_nodes;
  
   //number of time steps
  int num_dt = 4000;
  
  //number of processors
  int p;
  
  //local process array solution
  double *u1_local;
  
  //wall clock time
  double wtime;

  //Initialize MPI session
  MPI_Init(&argc, &argv);
  
  //returns the calling process's I.D.
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  //total number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  //tell user about the number processors used, the number of time steps and the total number of nodes used by the master process
  
  if (id == 0) 
  {
    timestamp();
    printf("Using %d processes.\n", p);
    printf("Using a total of %d points.\n", master_num_nodes);
    printf("Using %d time steps of size %g.\n", num_dt, dt);
  }

  //calculate wall clock time
  wtime = MPI_Wtime();

 //--Determine the number of nodes for local process

  i_master_min = (id*(master_num_nodes - 1))/p;
  i_master_max = ((id + 1)*(master_num_nodes - 1))/p;
  
  if (0 < id)
  {
    i_master_lo = i_master_lo - 1;
  }
  local_num_nodes = i_master_hi + 1 - i_master_lo;

  //update the local process nodes
  u1_local = update(id, p, master_num_nodes, local_num_nodes, num_dt, dt);
  
  //now collect the local values put them into master process/array
  collect(id, p, master_num_nodes, local_num_nodes, num_dt, dt, u1_local);

  //report wall clock time
  wtime = MPI_Wtime() - wtime;
  
  if (id == 0)
  {
    printf("\n");
    printf("Elapsed wallclock time was %g seconds\n", wtime);
  }

  //terminate the MPI session
  MPI_Finalize();

  //free the memory assigned to local process solution
  free(u1_local);

  //terminate the main function
  if (id == 0)
  {
    printf("\n");
    printf("Normal end of execution.\n");
    timestamp();
  }

  return 0;
}

//function to advance/move solution to the next time-step
double* update(int id, int p, int master_num_nodes, int local_num_nodes, int num_dt, double dt) 
{
  //Courant–Friedrichs–Lewy (CFL) condition
  double CFL; 
  
  //speed of wave
  double wave_speed;
  
  //spacial step between nodes
  double dx;
  
  //n represents the nodes in the time domain
  int n;
  
  //master array
  int i_master;
  
  //master array max value
  int i_master_hi;
  
  //master array min value
  int i_master_lo;
  
  //local array
  int i_local;
  
  //--local process arrays
  int i_local_max; //max value
  int i_local_min; //min value
  
  //used for transfer data between adjacent processes
  int ltor = 20; //left-to-right
  int rtol = 10; //right-to-left
  
  //MPI_Status represents the status of a reception operation
  MPI_Status status;
  
  //time time
  double t;
  
  //local solution arrays
  double *u0_local;
  double *u1_local;
  double *u2_local;
  
  //space point
  double x;
/*
  Determine the value of the Courant–Friedrichs–Lewy (CFL) condition:
  
  if 1 <= CFL, the solution is unstable; warn user.

  when 1 <= CFL, the time step is too large which means
  solutions won't be calculated at every node in the mesh

*/
  double wave_speed = 1.0;
  
  //spacial step
  dx = 1.0/(double)(master_num_nodes - 1);
  
  //now calculate the CFL condition
  CFL = wave_speed*dt/dx;

  //check for solution stability
  if (1.0 <= fabs(CFL))
  {
    if (id == 0)
    {
      fprintf(stderr,"\n");
      fprintf(stderr, "Courant-Friedrichs-Lewy Condition Update\n");
      fprintf(stderr, "Wave Speed = %g\n", c);
      fprintf(stderr, "dt = %g\n", dt);
      fprintf(stderr, "dx = %g\n", dx);
      fprintf(stderr, "CFL = %g\n", CFL);
      fprintf(stderr, "Computation will not be stable!\n");
    }
    
    //terminate MPI process
    MPI_Finalize();
    
    exit(1);
  }
/*
  
    to make use of our parallel strategy, domain decomposition, the master process
    must be divided into smaller processes.
    Each process stores about 1/p of the total + 2 extra slots
    
*/
  i_master_min = (id*(master_num_nodes - 1))/p;
  i_master_max = ((id + 1)*(master_num_nodes - 1))/p;
  
  if (0 < id)
  {
    i_master_min = i_master_min - 1;
  }

  i_local_min = 0;
  i_local_max = i_master_max - i_master_min;

  u0_local = (double*)malloc(local_num_nodes*sizeof(double));
  u1_local = (double*)malloc(local_num_nodes*sizeof(double));
  u2_local = (double*)malloc(local_num_nodes*sizeof(double));

  //starting at t = 0, we calculate the 1st spacial node value
  t = 0.0;
  
  //for the master process
  for (i_master = i_master_min; i_master <= i_master_max; i_master++) 
  {
    x = (double)(i_master)/(double)(master_num_nodes - 1);
    
    i_local = i_master - i_master_min;
    
    //local solution
    u1_local[i_local] = exact(x, t);
  }

  //for the local process
  for (i_local = i_local_min; i_local <= i_local_max; i_local++)
  {
    u0_local[i_local] = u1_local[i_local];
  }
  
  //time considerations. we have n nodes in the time domain
  //we now start iterating through all the time steps (4000)
  for (n = 1; n <= num_dt; i++)
  {
    t = dt*(double)n;

  //For the first time step, we need to use the initial derivative information.
    if (i == 1)
    {
      for (i_local = i_local_min + 1; i_local < i_local_max; i_local++) 
      {
        i_master = i_master_min + i_local;
        
        x = (double)(i_master)/(double)(master_num_nodes - 1);
        
        u2_local[i_local] = 
          
          + 0.5*CFL**2*u1_local[i_local-1]
          + (1.0 - CFL**2)*u1_local[i_local] 
          +  0.5*CFL**2*u1_local[i_local+1]
          +  dt*dudt(x, t);
      }
    }

  //now we can use the 2 previous solutions to calculate the next one

    else
    {
      for (i_local = i_local_min + 1; i_local < i_local_max; i_local++) 
      {
        u2_local[i_local] =
          
          + CFL**2*u1_local[i_local-1]
          + 2.0*(1.0 - CFL**2)*u1_local[i_local] 
          + CFL**2*u1_local[i_local+1]
          - u0_local[i_local];
      }
    }

    //Exchange data with "left-hand" neighbour/process
    if (0 < id) 
    {
      MPI_Send(&u2_local[i_local_min+1], 1, MPI_DOUBLE, id - 1, rtol, MPI_COMM_WORLD);
      MPI_Recv(&u2_local[i_local_min], 1, MPI_DOUBLE, id - 1, ltor, MPI_COMM_WORLD, &status );
    }
    else
    {
      x = 0.0;
      
      u2_local[i_local_min] = exact(x, t);
    }

    //Exchange data with "right-hand" neighbour/process
    if (id < p - 1) 
    {
      MPI_Send(&u2_local[i_local_max-1], 1, MPI_DOUBLE, id + 1, ltor, MPI_COMM_WORLD );
      MPI_Recv(&u2_local[i_local_max], 1, MPI_DOUBLE, id + 1, rtol, MPI_COMM_WORLD, &status);
    }
    else
    {
      x = 1.0;
      u2_local[i_local_max] = exact(x, t);
    }

    //Shift data for next time step.
    for (i_local = i_local_min; i_local <= i_local_max; i_local++)
    {
      u0_local[i_local] = u1_local[i_local];
      u1_local[i_local] = u2_local[i_local];
    }
  }

  //Free memory.
  free(u0_local);
  free(u2_local);

  return u1_local;
  
}

//send results from the local processes to the master process
void collect ( int id, int p, int master_num_nodes, int local_num_nodes, int num_dt, double dt, double u_local[]) 
{
  //buffer: region of a memory used to temporarily store data 
  //while it is being moved from 1 process to another 1
  int buffer[2];
  
  //--collect variables are used to store the master process index and values
  //collect1: contains the master index and values
  int collect1 = 10;
  //collect2: contains values for master process
  int collect2 = 20;
  
  //counter for iterating for workers
  int i;
  
  //master array
  int i_master;
  
  //--master array min and max values
  int i_master_max;
  int i_master_min;
  
  //local array
  int i_local;
  
  //local arrays min and max values
  int i_local_max;
  int i_local_min;
  
  //number of nodes for 2nd local process
  int local_num_nodes_2;
  
  //represents the status of a reception operation
  MPI_Status status;
  
  //time node value
  double t;
  
  //master solution array
  double *u_master;
  
  //spacial node value
  double x;

  //values for min and max values for master process arrays respectively
  i_master_min = (id*(master_num_nodes - 1))/p;
  i_master_max = ((id + 1)*(master_num_nodes - 1)/p;
                 
  if (0 < id)
  {
    i_master_min = i_master_min - 1;
  }

  //min and max values for local process arrays                
  i_local_lo = 0;
  i_local_max = i_master_max - i_master_min;

  //master process collects local results into the master solution array         
  if ( id == 0 )
  {

    //Create the global array.
    u_master = (double*)malloc(master_num_nodes*sizeof(double));

    //Copy the master's results into the global array.
    for (i_local = i_local_min; i_local <= i_local_max; i_local++)
    {
      i_master = i_master_min + i_local - i_local_min;
      u_master[i_master] = u_local[i_local];
    }

    //Contact each worker.
    for (i = 1; i < p; i++) 
    {

      //Message "collect1" contains the global index and number of values.
      MPI_Recv(buffer, 2, MPI_INT, i, collect1, MPI_COMM_WORLD, &status);
      
      //index of master process
      i_master_min = buffer[0];
      
      //nodes in 2nd process takes values from master process
      local_num_nodes_2 = buffer[1];

      if (i_master_lo < 0)
      {
        fprintf(stderr, "Illegal I_GLOBAL_LO = %d\n", i_global_min);
        exit ( 1 );
      }
      else if (master_num_nodes <= i_master_lo + n_local2 - 1)
      {
        fprintf(stderr, "  Illegal I_GLOBAL_LO + N_LOCAL2 = %d\n", i_master_min + local_num_nodes_2);
        exit(1);
      }

      //Message "collect2" contains the values.
      MPI_Recv(&u_master[i_master_min], local_num_nodes_2, MPI_DOUBLE, i, collect2, MPI_COMM_WORLD, &status);
    }
    
    //Print the results.
    t = dt*(double)num_dt;
    
    printf("\n");
    printf("    I      X     F(X)   Exact\n");
    printf("\n");
    
    for (i_master = 0; i_master < master_num_nodes; i_master++) 
    {
      x = (double)(i_master)/(double)(master_num_nodes - 1);
      
      printf("  %3d  %6.3f  %6.3f  %6.3f\n", i_master, x, u_master[i_master], exact(x, t));
    }

    free(u_master);
  }
                 
  //Workers send results to process 0.
  else
  {
    //Message "collect1" contains the global index and number of values.
    buffer[0] = i_master_min;
    buffer[1] = local_num_nodes;
    
    MPI_Send(buffer, 2, MPI_INT, 0, collect1, MPI_COMM_WORLD);

    //Message "collect2" contains the values.
    MPI_Send(u_local, local_num_nodes, MPI_DOUBLE, 0, collect2, MPI_COMM_WORLD);
  }

  return;
}

//function to evaluate exact solution
double exact(double x, double t)
{
  const double wave_speed = 1.0;
  
  #define pi 3.141592653589793
  
  //variable holds exact solution
  double value;

  value = sin(2.0*pi*(x - c*t));

  return value;
}

//for evaluating time derivative of the solution
double dudt(double x, double t)
{
  const double wave_speed = 1.0;
  
  #define pi 3.141592653589793
  
  double value;

  value = - 2.0*pi*c*cos(2.0*pi*(x - c*t));

  return value;
}

//for time considerations 
void timestamp()
{
#define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm* tm;
  time_t now;

  now = time(NULL);
  tm = localtime(&now);

  strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

  printf("%s\n", time_buffer);

  return;
  
#undef TIME_SIZE
}

