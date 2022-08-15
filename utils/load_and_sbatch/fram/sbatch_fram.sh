#!/bin/bash
#SBATCH --account=NN9561K 
#SBATCH --job-name=we0p02
#SBATCH --time=06-23:30:00
#SBATCH --nodes=32
#SBATCH --ntasks-per-node=32

set -o errexit
set -o nounset

module --quiet purge
module load intel/2018b

#cd /cluster/work/users/nicolos/hst/low_temp_1/hst_rk_we0p02/
cd /cluster/work/users/"your_name"/

mpirun -np 1024 ./flutas > my_output_file.txt 2>&1
