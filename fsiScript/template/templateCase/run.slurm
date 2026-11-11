#!/bin/bash
#SBATCH --job-name=dir_neo_1e-5_5e-4
#SBATCH --output=log.%x_%j.out
#SBATCH --error=log.%x_%j.err
#SBATCH --partition=long
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=168:00:00

# Good practice to navigate to the submission directory
cd $SLURM_SUBMIT_DIR

# Load OpenFOAM
source /usr/lib/openfoam/openfoam2412/etc/bashrc
#source /usr/lib/openfoam/openfoam2312/etc/bashrc
#source /usr/lib/openfoam/openfoam2212/etc/bashrc
#source /usr/lib/openfoam/openfoam2112/etc/bashrc

# Set number of threads for programs that use OpenMP parallelisation
# Note that OMP_NUM_THREADS defaults to 1 if not set
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export OMP_NUM_THREADS=10

# Execute the run script
./Allrun &> log.Allrun
