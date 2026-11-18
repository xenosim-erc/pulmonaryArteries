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

# Bind each MPI rank to one of the cores Slurm gave this job, in order. The
# cluster does not confine jobs to their cores (no cgroups), so unbound ranks
# drift onto cores busy with other jobs and run at half speed, holding back
# every other rank. CPU_IDs from scontrol are the OS core numbers here (one
# thread per core), which Open MPI's cpu list takes as given.
CPU_LIST=$(scontrol show job -d "$SLURM_JOB_ID" |
    sed -n 's/.*CPU_IDs=\([^ ]*\).*/\1/p' |
    head -n 1)
if [ -n "$CPU_LIST" ]
then
    export OMPI_MCA_hwloc_base_binding_policy=cpu-list:ordered
    export OMPI_MCA_hwloc_base_cpu_list="$CPU_LIST"
    echo "Binding MPI ranks to cores $CPU_LIST"
else
    echo "Could not read this job's cores from scontrol; MPI ranks are not bound" >&2
fi

# Report the cores each MPI rank may run on, sorted by rank: with the binding
# above, each rank should have one core, all different, together CPU_LIST
mpirun -np "$SLURM_NTASKS" python3 -c '
import os
rank = int(os.environ["OMPI_COMM_WORLD_RANK"])
cpus = sorted(os.sched_getaffinity(0))
print(f"Rank {rank}: CPUs {cpus}", flush=True)
' | sort -n -k2

# Execute the run script
./Allrun &> log.Allrun
