#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=4
#SBATCH --ntasks-per-node=4            # ranks per node, so 4 means each gpu gets a rank
#SBATCH --partition=debug              #compute_full_node, debug
#SBATCH --time=00:20:00                 # hours, minutes, seconds.  max time 24:00:00 

#SBATCH -o outjob_test.o%j       # Name of stdout output file
#SBATCH -e outjob_test.e%j       # Name of stderr error file
#SBATCH -J test

module load gcc cmake cuda openmpi
srun --mpi=pmix ./athena -i gamma1_50_beta1000_S100000_gpa53_33_nx6400_epsp0_1333.athinput
