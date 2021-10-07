#DIR="$(cd "$(dirname "$0")" && pwd)"
#export CXX=CC
#export CC=cc
export FC=mpiifort

#module load buildenv-intel/2018a-eb
#module load FFTW/3.3.6-nsc1
module restore system   # Restore loaded modules to the default
module load intel/2018b

make clean
make
