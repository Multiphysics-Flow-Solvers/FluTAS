## Synopsis

**FluTAS (Fluid Transport Accelerated Solver)** is an open-source code targeting multiphase fluid dynamics simulations. The key feature of **FluTAS** is the ability to efficiently run both on **many-CPUs** and **many-GPUs** in a single and unified framework. This framework is also designed to be modular so that it can be extended in a sustainable manner. So far, the code has been tested on several large-scale HPC clusters, both standard CPU-based and accelerated ones.

The overall effort was initiated within the research group of Prof. Luca Brandt at KTH in Stockholm (Sweden). The code structure and the Poisson solver are based on **CaNS** (https://github.com/p-costa/CaNS), and have been adapted to include the different modules developed by the various researchers and students at KTH Engineering Mechanics.

**References for the original code structure and Poisson solver (CPU and GPU)**

 * P. Costa. *A FFT-based finite-difference solver for massively-parallel direct numerical simulations of turbulent flows.* *Computers & Mathematics with Applications* 76: 1853--1862 (2018). [doi:10.1016/j.camwa.2018.07.034](https://doi.org/10.1016/j.camwa.2018.07.034) [[arXiv preprint]](https://arxiv.org/abs/1802.10323)

 * P. Costa, E. Phillips, L. Brandt & M. Fatica, *GPU acceleration of CaNS for massively-parallel direct numerical simulations of canonical fluid flows* *Computers & Mathematics with Applications* (2020). [doi.org/10.1016/j.camwa.2020.01.002](https://doi.org/10.1016/j.camwa.2020.01.002) [[arXiv preprint]](https://arxiv.org/abs/2001.05234)

**Reference for FluTAS and the acceleration of the multiphase code**

Crialesi-Esposito, M., Scapin, N., Demou, A. D., Rosti, M. E., Costa, P., Spiga, F., & Brandt, L. (2022). *FluTAS: A GPU-accelerated finite difference code for multiphase flows*. [[arXiv preprint]](https://arxiv.org/abs/2204.08834)

We recommend new users to have a look at this document and to carefully read the section **Compilation and usage**.

## News
 * not yet...to come!

## Code structure
To target different flow configurations, **FluTAS** has multiple `main.f90` each of one corresponding to a different application. Currently, we can target:

 * Single-phase solver, optionally with heat transfer;
 * Multiphase flows, using the VoF MTHINC method to capture the interface;
 * Multiphase flows with heat transfer, solving the energy equation (Boussinesq effects can be included in the liquid and gas phase).

In all the `main.f90`, different modules and subroutines are called and employed, depending on the specific problem under study. Although some apps are already present, more can be created by following the available templates. See the directory [`apps`](./src/apps/) for an overview of the existing applications.

Additional features are:

 * MPI parallelization in CPU
 * MPI+OpenACC+CUDA Fortran parallelization in GPU
 * FFTW guru interface used for computing multi-dimensional vectors of 1D transforms
 * The right type of transformation (Fourier, Cosine, Sine, etc) automatically determined from the input file
 * 2DECOMP&FFT routines used for performing global data transpositions and data I/O
 * Different physical configurations can be simulated just by selecting the appropriate application and adjusting the input files

Some examples of the flows that can be solved by FluTAS are:

 * periodic or spatially developing wall-bounded flows
 * tri-periodic domains, e.g., homogeneous isotropic turbulence
 * lid-driven cavity
 * Rayleigh–Bénard convection

The above configurations can be numerically reproduced both in single and multiphase flows, with or without heat transfer effects, in CPUs and in GPUs.
 
## Compilation and usage
Before compiling and using **FluTAS**, it is necessary to fulfill the prerequisites listed in [`REQ`](./getting_started/REQ.md). Once done, **FluTAS** can be compiled, run and the output can be visualized in Paraview. Details on how to perform all these tasks are provided as it follows:
 * For information about compilation and generic use, please refer to [`HOW_TO_USE`](./getting_started/HOW_TO_USE.md)
 * For information about the required input files, i.e., the files with extension .in, please refer to [`INFO_INPUT`](./getting_started/INFO_INPUT.md)
 * For information about how to visualize the output, please refer to [`INFO_VISU`](./getting_started/INFO_VISU.md)

## Notes
Please read the [`ACKNOWLEDGEMENTS`](./authorship/ACKNOWLEDGEMENTS.md) and [`LICENSE`](./authorship/LICENSE) files to have more information about the authors.

## Collaborate on this project
We highly value any feedback to improve the code and further extend it. In case of problems and bugs, please open an [issue](https://github.com/Multiphysics-Flow-Solvers/FluTAS/issues) to point out them and a [pull request](https://github.com/Multiphysics-Flow-Solvers/FluTAS/pulls) if you have a solution to fix them. Moreover, if you are interested in contributing to this effort, please get in contact with [Luca Brandt](mailto:luca@mech.kth.se).
