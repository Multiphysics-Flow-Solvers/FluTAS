## Prerequisites
**FluTAS** requires some external libraries for being compiled and run. In both CPU and GPU versions, the code parallelization relies on the library 2DECOMP, which is already included in the downloadable source code. Moreover, it requires an updated version of MPI and OpenMP. These are typically provided in standard supercomputers and can be easily installed in personal workstations and laptops.

### CPU version
For the use in ***CPU***, **FluTAS** requires the library FFTW to perform Fast Fourier transforms. FFTW should be downloaded, compiled separately and linked in the chosen `src/targets/*` file. To do so, we list here a series of steps:

1. Even before downloading/cloning **FluTAS**, create a separate directory in a location of your choice where to place the library. Any choice is fine but we strongly recommend creating an independent directory outside any numerical code you might have. On the command line, one can type `mkdir numerical_libraries` to create it, `cd numerical_libraries` to go inside it, and use the command `pwd` to know the exact path. Keep the information about the path since it will be important for steps 4 and 5;

2. Go to the page [download FFTW](http://www.fftw.org/download.html) and check which is the latest release of the library. While we are writing this note (August 2022), the latest release is `FFTW 3.3.10` and, therefore, we specify the procedure for this version. In case of subsequent releases, the procedure can be easily adjusted;
 
3. Go inside the new directory, i.e., `cd numerical_libraries`, and on the terminal type the following commands (or create a bash script to be run on the terminal)

   ~~~
   wget https://www.fftw.org/fftw-3.3.10.tar.gz
   tar xzf fftw-3.3.10.tar.gz
   cd fftw-3.3.10
   ~~~
 
4. Suppose your path saved at point 1 is `/home/numerical_libraries/fftw-3.3.10/`, on the terminal type the following commands (or create a bash script to be run on the terminal)

   ~~~
   ./configure FC=mpif90 CC=cc CXX=CC --prefix=/home/numerical_libraries/fftw-3.3.10/ --enable-threads --enable-openmp --enable-mpi
   make
   make install
   ~~~
 
5. Use the path where you have installed `fftw-3.3.10` (i.e., `/home/numerical_libraries/fftw-3.3.10/` in our case) to link **FluTAS** with the FFTW library in the chosen `src/targets/*`. You can simply do this by adding your path as an argument of `FFTW_HOME :=`. If there is an existing one, replace it with yours.

So far we have specified the compilation procedure, taking as an example the [`GNU compiler`](./../src/targets/target.generic-gnu). In case you plan to use alternative Fortran and/or C compilers, the compilation options `FC`, `CC` and `CXX` and the chosen `src/targets/*` should be changed accordingly.
 
### GPU version
For the use in ***GPU***, two NVIDIA libraries cuFFT and cuRAND from the CUDA toolkit are required and can be downloaded from the NVIDIA HPC SDK package (https://developer.nvidia.com/hpc-sdk). 


