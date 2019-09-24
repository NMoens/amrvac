#!/bin/sh
#SBATCH --job-name=atms3
#SBATCH --account=ivsusers
#SBATCH --time 1-00:00:00
#  (estimated run time in minutes)
# SBATCH --cpus-per-task=24
#SBATCH --tasks-per-node=24
#SBATCH -N 2
#  (default value, use this value if you want to execute a job on multiple nodes (openmpi))
#SBATCH --mem=2048
#  (memory in MB)
#SBATCH --partition=normal
#   (use the normal partition=default)
#SBATCH --output=/lhome/nicolasm/amrvac/mytests/OutflowAtmosphere3/nicolasm-stdout.log
#SBATCH --error=/lhome/nicolasm/amrvac/mytests/OutflowAtmosphere3/nicolasm-stderr.log

srun date

export AMRVAC_DIR=/lhome/nicolasm/amrvac
PATH="$AMRVAC_DIR:$AMRVAC_DIR/tools:./:$PATH"

module load mpi/openmpi-x86_64

cd /lhome/nicolasm/amrvac/mytests/OutflowAtmosphere3

make -j 8

mpiexec -N 32 amrvac -i usr_slurm.par 1>stdout.log 2>stderr.log

srun date
