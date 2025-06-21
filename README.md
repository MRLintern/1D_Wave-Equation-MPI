## 1D_Wave-Equation-MPI
* __1D Wave Equation__ Discretized using __Finite Differences__ and Solved via __Parallelizing__ the ___Jacobi Method___ with __MPI__.

## TODO
* Create a `Python` script which creates 2 sub-plots: `u(x)`, the wavefield and the error `u(x) - u_exact(x)`.

## Background: Parallel Programming

The parallel strategy used is ___Domain Decomposition___.
The problem (__global domain__) is decomposed into ___sub-domains___ (smaller processes). ___"Workers"___ in the sub-domains perform the calculations
and then communicate the results with the master (global domain). This [link](https://www.mcs.anl.gov/research/projects/mpi/tutorial/mpiexmpl/src2/io/C/main.html) provides a basic example behind the Master/Slave concept.

## Background: Wave Equation (A)

* The __Wave Equation__ takes the form:

        Utt = c^2*Uxx
* `Utt` is the 2nd time derivative of the `displacement field`.
* `Uxx` is the 2nd spacial derivative of the `displacement field`.

* The __1D version__ of the __Wave Equation__ is:

        d^2 u/dt^2 - c^2 * d^2 u/dx^2 = 0

        c is the speed of the wave, u is the displacement (field), t is time and x is the spatial component of the wave.

        space-interval: [x1, x2]
        time-interval: [t1, t2]

        Initial Conditions

        u(x,t1) = u_t1(x); u(x,0) = g(x,t=0) = sin(2*pi*(x-c*t))
        u(x,t1) = ut_t1(x); dudt(x,0) = h(x,t=0) = -2*pi*c*cos(2*pi*(x-c*t))
        
* Boundary (Dirichlet) conditions: 

        u(x1,t) = u_x1(t); u(0,t) = u0(t) = sin(2*pi*(0-c*t))
        u(x2,t) = u_x2(t); u(1,t) = u1(t) = sin(2*pi*(1-c*t))

* Discretized version of the wave equation:

        uxx = (u(x+dx,t) - 2u(x,t) + u(x-dx,t))/dx^2
        utt = (u(x,t+dt) - 2u(x,t) + u(x,t-dt))/dt^2

* After some algebra and simplification, we end up with the final finite difference equation:

        u(i,n+1) = -u(i,n-1) + 2u(i,n) + CFL^2(u(i+1,n) - 2u(i,n) + u(i-1,n)),

where `n` represents the nodes in the time direction and `i` represents the nodes in the spacial direction.

## Background: Wave Equation (B)

* In the last equation, the term `CFL` was introduced. This represents the `Courant–Friedrichs–Lewy condition`.
* It takes the form:

        CFL = c*dt/dx,

where `c` is the wave speed, `dx` the spacial step and `dt` the time step.
* This condition ensures that every node in the discretized mesh is calculated.
* The value of `CFL` needs to be in range: 

        0.5 < CFL < 1

* If we have:

        CFL > 1,

then nodes will be missed and the solution will not be smooth and becomes unstable leading to errors. 



## Requirements

* Operating Systems: `Ubuntu 20.04`.
* Compiler: `mpicc`.
* `MPICC` so `MPI` can be used.
* Text Editor. Any can be used. E.g. `Visual Studio Code`, `Vim`, `Emacs`, `gedit` etc.

## Installing MPI

`MPI`can be downloaded [here](https://www.mpich.org/) or at the command line (see below).

* `$ sudo apt-get install mpich`

## Running the Application

* The number of required processes to be used can be specified by `-np` (see below), where 4 have been used. The compiler optimisation flag in the `Makefile` has been set to `-O3`. This can be changed of course.
* You can run `make test` to run the `Unit Test Script` prior to building the final software.
* This will generate a results file (`computed vs. analytic data`) called `wave_test_output.txt` within the `results` directory; the `results` directory is built by `Make`.

* `$ git clone https://github.com/MRLintern/1D_Wave-Equation-MPI.git`
* `$ make`
* My machine has `4 CPUs`, so change accordingly.
* `$ mpirun -np 4 ./main`
* To 'clean up' the application run: `$ make clean`. This will get rid of `bin`, `results/results.txt` and  `wave_test_output.txt`.
* Note: the executable `main` in the parent directory will still be present. This is due to `main.c` being compiled in the `tests` (for `unit testing`) and parent directories.


## Results

* The results of the calculations can be found in the `results` folder.
* 1, 2 and 4 processors were used.
* Note: this is an __example folder__. When you `make` the application, a folder called `results` is generated with the output data.

