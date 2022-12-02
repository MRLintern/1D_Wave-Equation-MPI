## Note:

Software is currently being redeveloped. Once its all sorted, it'll appear here.

## 1D_Wave-Equation-MPI
1D Wave Equation Discretized using Finite Differences and Solved via Parallelizing the Jacobi Method with MPI

## Requirements

* Operating Systems: `Ubuntu 20.04`.
* Compiler: `gcc 9.4.0`.
* `mpich` so `MPI` can be used.
* Text Editor. E.g. `Visual Studio Code`, `Vim`, `Emacs`, `gedit` etc.

## Installing MPI

[MPI](https://www.mpich.org/) can either be downloaded via the link provided or at the command-line.

* `$ sudo apt-get install mpich`

## Running the Application

* `$ git clone https://github.com/MRLintern/1D_Wave-Equation-MPI.git`
* `$ chmod +x wave.sh`
* `$ ./main_w`

## Results

The results of the calculations are printed to `results.txt`. 

## Additional

Vary the number of processors used and the different compiler optimisation flags. 
A plot with the different speed-up values using different numbers of processors and compiler optimisation flags is included.
