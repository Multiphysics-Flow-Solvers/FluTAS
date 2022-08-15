#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --ntasks-per-socket=2
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:4
#######SBATCH --mem=230000MB
#SBATCH --time 00:05:00
#SBATCH -A pra22_5578
#SBATCH -p    m100_usr_prod
#######SBATCH --qos=m100_qos_bprod
#SBATCH --job-name=FluTAS_job_test_1
#SBATCH --error=log.%j.err
#SBATCH --output=log.%j.out

#SBATCH --mail-type=ALL
#######SBATCH --mail-user=marco.crialesiesposito@gmail.com
#SBATCH --mail-user=demou@mech.kth.se

module purge
module load autoload profile/global
module load  fftw/3.3.8--spectrum_mpi--10.3.1--binary
module load hpc-sdk/2020--binary

#export OMP_NUM_THREADS=1
#export OMP_PROC_BIND=true
#export NO_STOP_MESSAGE=yes
#export PGI_ACC_TIME=1
#export PGI_ACC_NOTIFY=2
#export CUDA_VISIBLE_DEVICES="$SLURM_LOCALID % $SLURM_GPUS_PER_NODE"
LD_LIBRARY_PATH="/cineca/prod/opt/compilers/cuda/10.0/none/extras/CUPTI/lib64/:$LD_LIBRARY_PATH"

mpirun -np 16 --map-by socket:PE=8 --rank-by core --report-bindings ./flutas

