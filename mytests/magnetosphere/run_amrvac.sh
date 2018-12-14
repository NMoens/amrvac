#!/bin/sh
#SBATCH --job-name=msph
#SBATCH --account=ivsusers
#SBATCH --time 14-00:00:00
#  (estimated run time in minutes)
# SBATCH --cpus-per-task=24
#SBATCH --tasks-per-node=24
#SBATCH -N 2
#  (default value, use this value if you want to execute a job on multiple nodes (openmpi))
#SBATCH --mem=2048
#  (memory in MB)
#SBATCH --partition=long
#   (use the normal partition=default)
#SBATCH --output=/STER/nicolasm/amrvac2/amrvac/mytests/magnetosphere/nicolasm-stdout.log
#SBATCH --error=/STER/nicolasm/amrvac2/amrvac/mytests/magnetosphere/nicolasm-stderr.log

srun date

export AMRVAC_DIR=$HOME/../../STER/nicolasm/amrvac2/amrvac 
PATH="$AMRVAC_DIR:$AMRVAC_DIR/tools:./:$PATH"

module load mpi/openmpi-x86_64

cd /STER/nicolasm/amrvac2/amrvac/mytests/magnetosphere

mpiexec amrvac -i usr_slurm.par 1>stdout.log 2>stderr.log

srun date
