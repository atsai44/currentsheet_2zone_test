#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-node=1           # ranks per node, so 4 means each gpu gets a rank
#SBATCH --partition=compute
#SBATCH --time=04:00:00                 # hours, minutes, seconds

module load python
srun python 2D_plots.py
