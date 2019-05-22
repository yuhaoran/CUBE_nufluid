#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=80
#SBATCH --time=1:00:00
#SBATCH --job-name test_main
#SBATCH --output=output_main.txt

cd $SLURM_SUBMIT_DIR
source ../utilities/module_load_niagara_intel.sh
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

env > env.txt
ulimit -a > ulimit.txt
srun -N 1  ./main.x
#mpirun -bind-to node -np 1 ./main.x 
