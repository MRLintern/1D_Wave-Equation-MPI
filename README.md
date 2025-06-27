## 1D_Wave-Equation-MPI
---
* This project models the __1D wave equation__ using __MPI__ (__Message Passing Interface__) for parallel computation.
* The main goal is to demonstrate how __distributed-memory parallelism__ can be used to solve __partial differential equations (PDEs)__ like the wave equation efficiently.
## Update:
---
* The software is currently being refactored into Modern and modular form.
## Background: Parallel Programming
---
* The parallel strategy used is ___Domain Decomposition___.
* The problem (__global domain__) is decomposed into ___sub-domains___ (smaller processes). ___"Workers"___ in the sub-domains perform the calculations
and then communicate the results with the master (global domain). This [link](https://www.mcs.anl.gov/research/projects/mpi/tutorial/mpiexmpl/src2/io/C/main.html) provides a basic example behind the Master/Slave concept.
* Each process updates its section of the __domain independently__ but must __exchange__ __"ghost" cell data__ to compute __spatial derivatives__.

#### Why Use Domain Decomposition?

* Reduces the memory usage per process.
* Enables __concurrent computations__.
* Minimises idle time by __overlapping computatiion__ and __communication__.


## The Wave Equation
---
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
        
* Boundary (__Dirichlet__) conditions: 

        u(x1,t) = u_x1(t); u(0,t) = u0(t) = sin(2*pi*(0-c*t))
        u(x2,t) = u_x2(t); u(1,t) = u1(t) = sin(2*pi*(1-c*t))

* Discretized version of the wave equation:

        uxx = (u(x+dx,t) - 2u(x,t) + u(x-dx,t))/dx^2
        utt = (u(x,t+dt) - 2u(x,t) + u(x,t-dt))/dt^2

* After some algebra and simplification, we end up with the final __finite difference equation__:

        u(i,n+1) = -u(i,n-1) + 2u(i,n) + CFL^2(u(i+1,n) - 2u(i,n) + u(i-1,n)),

where `n` represents the nodes in the time direction and `i` represents the nodes in the spacial direction.

## The Courant–Friedrichs–Lewy Condition
---
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
---
* Operating Systems: `Ubuntu 20.04`.
* Compiler: `mpicc`.
* `MPICC` so `MPI` can be used.
* Text Editor. Any can be used. E.g. `Visual Studio Code`, `Vim`, `Emacs`, `gedit` etc.

## Installing MPI
---
`MPI`can be downloaded [here](https://www.mpich.org/) or at the command line using Ubuntu's package manager:

* `$ sudo apt-get install mpich`

## Running the Application
---
### Unit Testing
* A smaller version of the software can be found in the `tests` directory for unit testing purposes.
* No results are produced, only the functions `update()` and `collect()`, which are tested.
* The software was tested on `2 processors`.
#### Instructions for Unit Testing:
* `$ make test`
* `$ mpirun -np 2 ./test`
* Once happy, clean up the directory. I.e. get rid of the binary file and the `.o` files.
* `$ make clean`
* TODO: New instructions to be added once the new software has been developed.

## Results
---
* TODO: new results to be generate; a `Python` script which plots the `computed u(x)` and the `u_exact(x)` values for processor numbers; `1,2,3 & 4`. 

