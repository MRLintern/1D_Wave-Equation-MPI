## 1D_Wave-Equation-MPI
1D Wave Equation Discretized using Finite Differences and Solved via Parallelizing the Jacobi Method with MPI

## Background: Parallel Programming

The parallel strategy used is Domain Decomposition.
The problem (global domain) is decomposed into smaller domains (processes). "Workers" in the subdomains perform the calculations
and then communicate the results with the master (global domain). This [link](https://www.mcs.anl.gov/research/projects/mpi/tutorial/mpiexmpl/src2/io/C/main.html) provides a basic example behind the Master/Slave concept.

## Background: Wave Equation (A)

The Wave Equation takes the form:

        Utt = c^2*Uxx

1D version of the Wave Equation:

        d^2 u/dt^2 - c^2 * d^2 u/dx^2 = 0

        c is the speed of the wave, U is the displacement (field), t is time and x is the spatial component of the wave.

        space-interval: [x1, x2]
        time-interval: [t1, t2]

        Initial Conditions

        U(x,t1) = U_t1(x); u(x,0) = g(x,t=0) = sin(2*pi*(x-c*t))
        U(x,t1) = Ut_t1(x); dudt(x,0) = h(x,t=0) = -2*pi*c*cos(2*pi*(x-c*t))

        
Boundary (Dirichlet) conditions: 

        U(x1,t) = U_x1(t); u(0,t) = u0(t) = sin(2*pi*(0-c*t))
        U(x2,t) = U_x2(t); u(1,t) = u1(t) = sin(2*pi*(1-c*t))

Discretized version of the wave equation:

        Uxx = (U(x+dx,t) - 2U(x,t) + U(x-dx,t))/dx^2
        Utt = (U(x,t+dt) - 2U(x,t) + U(x,t-dt))/dt^2

After some algebra and simplification, we end up with the final finite difference equation:

        u(i,n+1) = -u(i,n-1) + 2u(i,n) + alpha^2(u(i+1,n) - 2u(i,n) + u(i-1,n))

where `n` represents the nodes in the time direction and `i` represents the nodes in the spacial direction.

## Background: Wave Equation (B)

In the last equation, the term `alpha` was introduced. This represents the `Courant–Friedrichs–Lewy condition`.
It takes the form:

        alpha = c*dt/dx,

where `c` is the wave speed, `dx` the spacial step and `dt` the time step.
This condition ensures that every node in the discretized mesh is calculated.
The value of `alpha` needs to be in range: 

        0.7 < alpha < 1

If we have:

        alpha => 1,

then nodes will be missed and the solution will not be smooth and becomes unstable leading to errors. 



## Requirements

* Operating Systems: `Ubuntu 20.04`.
* Compiler: `gcc 9.4.0`.
* `mpich` so `MPI` can be used.
* Text Editor. Any can be used. E.g. `Visual Studio Code`, `Vim`, `Emacs`, `gedit` etc.

## Installing MPI

`MPI`can be downloaded [here](https://www.mpich.org/) or at the command line (see below).

* `$ sudo apt-get install mpich`

## Running the Application

* `$ git clone https://github.com/MRLintern/1D_Wave-Equation-MPI.git`
* `$ chmod +x wave.sh`
* `$ ./main_w`
* or `$ bash main_w`

## Results

The results of the calculations are printed to `results.txt`. 

## Additional

Vary the number of processors used and the different compiler optimisation flags. 
A plot with the different speed-up values using different numbers of processors and compiler optimisation flags is included.
