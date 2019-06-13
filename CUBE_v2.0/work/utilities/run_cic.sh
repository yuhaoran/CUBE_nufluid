#!/bin/bash
#SBATCH --job-name=caf
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=8GB

source ../utilities/module_load_prince.sh
make clean
make


srun ./cicpower.x > log_cicpower
