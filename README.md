## Synopsis

**FluTAS (Fluid Transport Accelerated Solver)** is a modular, multiphysics code for multiphase fluid dynamics simulations. The code is written in a modular way, to easily accomodate different physics into the system. One of the main purposes of the project is to provide an efficient framework able to run both on many-CPUs (MPI) and many-GPUs (MPI+OpenACC+CUDA Fortran). The code has been tested on several large-scale HPC clusters.

This effort initiated within the research group of Prof. Luca Brandt at KTH. The code structure and the Poisson solver are based on **CaNS** (https://github.com/p-costa/CaNS), and have been adapted to include the different modules developed by the various researchers and students at KTH Engineering Mechanics.

**References**

**Method and single-phase solver**:  
P. Costa. *A FFT-based finite-difference solver for massively-parallel direct numerical simulations of turbulent flows.* *Computers & Mathematics with Applications* 76: 1853--1862 (2018). [doi:10.1016/j.camwa.2018.07.034](https://doi.org/10.1016/j.camwa.2018.07.034) [[arXiv preprint]](https://arxiv.org/abs/1802.10323)

**GPU extension**:  
P. Costa, E. Phillips, L. Brandt & M. Fatica, *GPU acceleration of CaNS for massively-parallel direct numerical simulations of canonical fluid flows* *Computers & Mathematics with Applications* (2020). [doi.org/10.1016/j.camwa.2020.01.002](https://doi.org/10.1016/j.camwa.2020.01.002) [[arXiv preprint]](https://arxiv.org/abs/2001.05234)

## News
 * not yet...to come!

## Features
Current physics modules implemented:
 * Base solver for the incompressible two-fluid Navier-Stokes equations;
 * Multiphase flows, using the VoF MTHINC method;
 * Heat Transfer, solving the energy equation (Boussinesq effects can be included in the gas phase). 

Some features are:

 * MPI parallelization
 * FFTW guru interface used for computing multi-dimensional vectors of 1D transforms
 * The right type of transformation (Fourier, Cosine, Sine, etc) automatically determined from the input file
 * 2DECOMP&FFT routines used for performing global data transpositions and data I/O
 * A different canonical flow can be simulated just by changing the input files

Some examples of the flows that can be solved by FluTAS are:

 * periodic or spatially developing channel
 * periodic or spatially developing square duct
 * tri-periodic domain
 * lid-driven cavity
 * homogeneous and isotropic turbulence
 * Rayleigh–Bénard convection
 
## Usage
for more details, please refer to [`HOW_TO_USE.md`](./HOW_TO_USE.md)

### Compilation
The code should be compiled in `src/`. The prerequisites are the following:

 * MPI
 * FFTW3
 * OpenMP (optional)

and for the GPU version:

 * NVIDIA Fortran compiler (https://developer.nvidia.com/hpc-sdk)
 * cuFFT from the CUDA toolkit

Further technical informations are provided in [`HOW_TO_USE.md`](./HOW_TO_USE.md).

### Running the code
Run the executable with `mpirun` with a number of tasks and shared threads complying to what has been set in the input file `dns.in`. Data will be written by default in a folder named `data/`, which must be located where the executable is run.

### Visualizing field data
See [`INFO_VISU.md`](./INFO_VISU.md).



## Notes
We appreciate any feedback to improve the code. Also, feel free to send case files pertaining to flows not listed in the examples folder.
Please read the [`ACKNOWLEDGEMENTS.md`](./ACKNOWLEDGEMENTS.md) and [`LICENSE`](./LICENSE) files.

## Collaborate on this project
If you are interested in contributing to this effort, please get in contact with [Luca Brandt](mailto:luca@mech.kth.se).
