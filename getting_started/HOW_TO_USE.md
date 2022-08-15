### First steps with **FluTAS**
We list some basic steps for beginner users of the code:

 1. In the first step, ensure that the prerequisites listed [`here`](./REQ.md) are fulfilled on the machine where you are going to run the code (either laptop, workstation or supercomputer);
 2. Clone or fork the code. To fork it, you can follow the explanation provided [`here`](https://docs.github.com/en/get-started/quickstart/fork-a-repo). For cloning it, you can simply type on the command line:
    ~~~
    git clone https://github.com/Multiphysics-Flow-Solvers/FluTAS.git
    ~~~
 3. Go inside the directory `src/` of **FluTAS**, i.e., `cd src/`. Now you are ready to compile it.
 
### Compiling **FluTAS**
The compilation step depends on which [`application`](./../src/apps) is chosen and on which [`architecture`](./../src/targets) the user is working. Currently, the code supports the following applications:

 1. `single_phase`: single-phase, incompressible, optionally with heat transfer effect.\
    To create the executable for this application, type on the command line:
    ~~~
    make clean APP=single_phase && make ARCH=generic-gnu APP=single_phase DO_DBG=0 -j4
    ~~~
 2. `two_phase_inc_isot`: two-phase, incompressible and isothermal flow.\
    To create the executable for this application, type on the command line:
    ~~~
    make clean APP=two_phase_inc_isot && make ARCH=generic-gnu APP=two_phase_inc_isot DO_DBG=0 -j4
    ~~~
 3. `two_phase_ht`: two-phase, incompressible flow with heat transfer.\
    To create the executable for this application, type on the command line:
    ~~~
    make clean APP=two_phase_ht && make ARCH=generic-gnu APP=two_phase_ht DO_DBG=0 -j4
    ~~~

Note that in the above example, we have compiled the code in parallel (the option `-j4` means using 4 processors). Moreover, we employ the GNU compiler to build code. Other compilation options can be found in the [`targets`](./../src/targets) folder, where INTEL, NVF and CRAY targets are available and can be simply linked changing the argument of `ARCH` above. For code development and extension, we recommend compiling **FluTAS** in debugging mode, i.e., setting `DO_DBG=1`. Furthermore, to ensure a correct compilation of the code, whether an application is chosen, a default `basic` choice is provided in apps/ which corresponds to the application `two_phase_inc_isot`. Note that in case of a misspelling of the chosen application, the default `basic` choice will be compiled.

If the compilation has been successfully performed, the next step is to fill the input files depending on the desired applications. In the [`examples`](./../examples) folder, there are examples of input files for several canonical flows, and the users can copy these files to the same directory where the executable is placed. We also recommend seeing [`INFO_INPUT`](./INFO_INPUT.md) for a detailed description of how the input files should be filled.

### Running **FluTAS**
Once the code is successfully compiled and the input files properly filled, **FluTAS** can be run. We recommend the following steps:
 1. Create a run folder outside the `src` directory, i.e., `mkdir run`. In your workstation, any choice is fine while on clusters we recommend cloning the code and compile it in your home directory (typically backed-up but with limited storage) and placing the executables (both `flutas` and `flutas.<YOUR_APPLICATION>`) together with the input files in your project folder (typically not backed-up but with more storage)
 2. On the command line (or in submission script), type (or write) `mpirun -n NP flutas` where `NP` is the number of processors or GPUs prescribed in dns.in, i.e., the product between the components of `dims_in` (set in the input file `dns.in`, see [`INFO_INPUT`](./INFO_INPUT.md)). Note that in case of different compilation choices (e.g., INTEL, CRAY), the command `mpirun` requires to be adjusted depending on the environment you are working and the compilation options you have chosen.

### Data visualization and Input/Output directives
The data generated during the execution of the program are stored in the directory `data`, automatically created by **FluTAS** in the same location where the executable is placed. Before running the simulation, it is important to decide what to print. To this end, the files `out1d.h90`, `out2d.h90` and `out3d.h90` in `src/apps/<YOUR_APPLICATION>/postp.<YOUR_APPLICATION>/` set which data are written in 1-, 2- and 3-dimensional output files, respectively. Replace `<YOUR_APPLICATION>` with the name of your application. Note that the code should be recompiled after editing `out*d.h90` files. The generated binary files can be read and visualized following the instruction reported [`here`](./INFO_VISU.md).

### Pre-processor flags
The code is compiled using different preprocessor flags which also control the employed modules and subroutines. If shared among different applications, these flags are in the Makefile. If they are specific to a certain application, they are typically placed in `src/apps/<YOUR_APPLICATION>/apps.<YOUR_APPLICATION>/`. As usual, remember to replace `<YOUR_APPLICATION>` with the name of your application. 

Currently, the following preprocessor flags can be used:
 * `-D_USE_CUDA`:                  enabled by setting `USE_CUDA=1`. It uses the variable `CUDA_LIB` to locate CUDA libraries (please adapt it to your system). Also, remember to adapt the compilation flags to your system (currently `-gpu=cc70,cuda11.0`).
 * `-D_USE_NVTX`:                  enables the NVTX profiling. Please modify the variable `LIB_NVTX` to your system.
 * `-D_CONSTANT_COEFFS_POISSON`:   enabled by setting `CONSTANT_COEFFS_POISSON=1`. **Please keep it to 1**. Future developments will include other variable coefficient Poisson solvers.
 * `-D_USE_VOF`:                   enables the VoF module by setting `USE_VOF=1`. _Reminder: it forces the usage of vof.in input file_ (see below). 
 * `-D_VOF_DBG`:                   Skips the flow solution and imposes a pure advection at constant velocity for VoF debugging purposes. 
 * `-D_HEAT_TRANSFER`:             enables the computation of the energy equation by setting `HEAT_TRANSFER=1`. _Reminder: it forces the usage of heat_transfer.in input file_ (see below). 
 * `-D_INIT_MONTECARLO`:           enabled with `INIT_MONTECARLO=1`, solves hyperbolic tangent problem for the consistent initialization of MTHINC VOF. **not tested in GPU**
 * `-D_DO_POSTPROC`:               enabled with `DO_POSTPROC=1`, allows for the computation of statistics and data extraction on-the-fly.  _Reminder: it forces the usage of post.in input file_ (see below). 
 * `-D_TURB_FORCING`:              enabled with `TURB_FORCING=1`, enables the external turbulence forcing, either with ABC or TGV methods.  _Reminder: it forces the usage of forcing.in input file_ (see below). 
 * `-D_TIMING`:                    enables timestep timing
 * `-D_TWOD`:                      enables the computation of 2D cases. The X direction is not computed and it should be set with periodic boundaries and 2 grid points.
 * `-D_BOUSSINESQ`:                enables the solution of the heat transfer equation using the Boussinesq approximation in the gas phase

### Further extensions
We are happy if researchers choose **FluTAS** as a base solver for their applications. For future and further developments, we provide the following recommendations:
 1. Clone or fork the repository (private or public, depending on the granted access you have);
 2. Work on a separate branch from the master or main one. To do so, type on the command line `git checkout -b "<NAME_OF_YOUR_BRANCH>"`. Choose a name for your meaningful name for your branch and work on it;
 3. You can use the already available applications or create new ones. In this case, we recommend following the already available [`templates`](./../src/apps);
 4. If you plan to create new applications:
     * create a new directory with the name of your app in [`apps`](./../src/apps) folder, i.e., `cp -r "<NAME_OF_AN_EXISTING_APP>" "<NAME_OF_YOUR_NEW_APP>"`. As an existing app, choose the one closest to your application and needs;
     * adjust the name of the existing files and folders inside the new directory. You can type on the terminal, `mv main__"<NAME_OF_OLD_APP>".f90 main__"<NAME_OF_YOUR_NEW_APP>".f90`, `mv post."<NAME_OF_OLD_APP>" post."<NAME_OF_YOUR_NEW_APP>"` and `mv app."<NAME_OF_OLD_APP>" app."<NAME_OF_YOUR_NEW_APP>"`;
     * open the file `app."<NAME_OF_YOUR_NEW_APP>"`, e.g., `vi app."<NAME_OF_YOUR_NEW_APP>"` and modify the name of the existing `main__"<NAME_OF_OLD_APP>".f90` with `main__"<NAME_OF_YOUR_NEW_APP>".f90`. In the same file, there is also a list of preprocessor flags specific to this app. Adjust them or add new ones depending on your needs;
     * modify the descriptions and comments inside `param.f90` and `main__"<NAME_OF_YOUR_NEW_APP>".f90`
     * try to create an executable for the new application to test if everything has been done properly
       ~~~
       make clean APP=<NAME_OF_YOUR_NEW_APP> && make ARCH=generic-gnu APP=<NAME_OF_YOUR_NEW_APP> DO_DBG=0 -j4
       ~~~
 5. Once the new application is created, the user can start to modify the source code. Many code styles are possible, but a modular and sustainable programming practice is strongly encouraged, i.e.,
     * Keep the subroutines as pure as possible with most of the variables declared with a defined intent, i.e., `in`, `out` or `inout`. This **input/output** approach must be respected for the "large" arrays (e.g., velocity, pressure, temperature fields, etc.)
     * Do not define global arrays visible to all the subroutines, but define all the arrays locally and only in the subroutines where they are needed. This is beneficial also for the GPU porting and to ensure efficient use of the available memory;
     * Add comments to the code and respect a minimum indentation for readability;
     * Do not make abuse of preprocessor macros, sometimes duplicate the subroutines rather than use too many preprocessor macros;
     * For each new development and new application, create different benchmarks of increasing complexities that serve as examples for the new application and to test, debug and continuously integrate the code. In the [`examples`](./../examples) directory, create a new folder for the new app and place the examples there. Ideally, just changing the input files and after having generated the executable for that application, the user is able to run a certain example without modifying the source code.

