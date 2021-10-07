### Code preparation 
Different functions or parts of the code can be used, depending on the specific problem under study. Therefore, users are encouraged to adapt the desired `main*.f90` file from the ones available in the folder `src/apps/`. Although some apps are already present, more  can be created by following the available templates.

### Code compilation
The code is compiled using the file `src/Makefile` appropriately configured. An example of compilation command is `make -j4 ARCH=generic APP=generic`. For further details on the `ARCH` and `APP` variables, please refer to the following sections.

### Architecture target
In the folder `src/targets/` several templates for architecture-depending environments can be found. These files includes examples for CPU and GPU clusters, as well as various compilers. The specific location where external libraries are installed can be set in the files `src/tagets/target*`, alongside the desired compilation flags.

### Application selection 
Different applications can be chosen to better reproduce the physical problem under study. The available applications can be found in `src/apps/`, where each one presents a `main*.f90` and a `app.*` file containing the preprocessor flags required to successfully compile and run the code.
Currently, the code supports the following applications:
	* `basic`: single phase, isothermal flow;
	* `two\_phase\_inc\_isot`: two-phase, incompressible and isothermal flow;
	* `two\_phase\_ht`: two-phase, incompressible flow with heat transfer.

### Pre-processor flags
While several `main.f90` can be produced, one `Makefile` is generally used for compilation. It is advisable to duplicate the `Makefile` within each folder corresponding to each configuration, to avoid inconsistent compilation.
Currently, the following preprocessor flags can be used:
 * `-D_USE_CUDA`:                  enabled by setting `USE_CUDA=1`. It uses the variable `CUDA_LIB` to locate CUDA libraries (please adapt it to your system). Also, remember to adapt the compilation flags to your system (currently `-gpu=cc70,cuda11.0`).
 * `-D_USE_NVTX`:                  enables the NVTX profiling. Please modify the variable `LIB_NVTX` to your system.
 * `-D_CONSTANT_COEFFS_POISSON`:   enabled by setting `CONSTANT_COEFFS_POISSON=1`. **Please keep it to 1**. Future developments will include other variable coefficient Poisson solvers.
 * `-D_USE_VOF`:                   enables the VoF module by setting `USE_VOF=1`. _Reminder: it forces the usage of vof.in input file_ (see below). 
 * `-D_VOF_DBG`:                   Skips the flow solution and imposes a pure advection at constant velocity for VoF debugging purposes. 
 * `-D_HEAT_TRANSFER`:             enables the computation of the energy equation by setting `HEAT_TRANSFER=1`. _Reminder: it forces the usage of heat_transfer.in input file_ (see below). 
 * `-D_INIT_MONTECARLO`:           enabled with `INIT_MONTECARLO=1`, solves hyperbolic tangent problem for tthe consistent initialization of MTHINC VOF. **not tested in GPU**
 * `-D_DO_POSTPROC`:               enabled with `DO_POSTPROC=1`, allows for the computation of statistics and data extraction on-the-fly.  _Reminder: it forces the usage of post.in input file_ (see below). 
 * `-D_TURB_FORCING`:              enabled with `TURB_FORCING=1`, enables the external turbulence forcing, either with ABC or TGV methods.  _Reminder: it forces the usage of forcing.in input file_ (see below). 
 * `-D_TIMING`:                    enables timestep timing
 * `-D_TWOD`:                      enables the computation of 2D cases. The X direction is not computed and it should be set with periodic boundaries and 2 grid points.
 * `-D_MOM_QUICK`:                 enables the QUICK scheme for the discretization of the advection term
 * `-D_BOUSSINESQ`:                enables the solution of the heat transfer equation using the Boussinesq approximation in the gas phase

### External libraries
Few external libraries are required for the code to compile and run. In both CPU and GPU versions, the code parallelization relies on the library 2DECOMP, which is already included in the downloadable code snapshot. For the CPU version of the code, FFTW should be compiled separately and linked in the respective `src/tagets/target*` file. 

The GPU version of the code uses the NVIDIA libraries already available in the NVIDIA HPC SDK (https://developer.nvidia.com/hpc-sdk) package, namely CUFFT and CURAND.

### Input files
The input files templates can be found in `utils/templates/` folder. In each of them, the corresponding lines are commented with the equivalent of the parameter name in the code (see file `param.f90`).

The input file dns.in sets the physical and computational parameters. In the `examples/` folder are examples of input files for several canonical flows. See `src/INFO\_INPUT.md` for a detailed description of the input file.

Files `out1d.h90`, `out2d.h90` and `out3d.h90` in `src/` set which data are written in 1-, 2- and 3-dimensional output files, respectively. The code should be recompiled after editing `out*d.h90` files.


