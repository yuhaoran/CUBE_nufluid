#!/bin/bash
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80
#SBATCH --time=0:25:00
#SBATCH --job-name test_main2
#SBATCH --output=output_main.txt

cd $SLURM_SUBMIT_DIR
source ../utilities/module_load_niagara_intel.sh
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export FOR_COARRAY_NUM_IMAGES=8
export OMP_NUM_THREADS=80

#cd ../utilities/
#./ic.x

cd ../main/
./main.x

#cd ../utilities/
#srun -N 1 ./cicpower.x

#mpirun -bind-to node -np 1 ./main.x 
