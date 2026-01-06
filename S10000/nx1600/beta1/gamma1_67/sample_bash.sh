#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=4              # 4 gpus per node
#SBATCH --ntasks-per-node=4            # ranks per node, so 4 means each gpu gets a rank
#SBATCH --partition=debug              #compute_full_node
#SBATCH --time=0:30:00                 # hours, minutes, seconds.  max time 24:00:00 , think its like 2 hours for 1 node half hour for 2

#SBATCH -o outjob_test.o%j       # Name of stdout output file
#SBATCH -e outjob_test.e%j       # Name of stderr error file
#SBATCH -J test

module load gcc cmake cuda openmpi
srun --mpi=pmix ./athena -i gamma1_67_beta1_S10000_gpa13_33_nx1600_epsp0_0000.athinput
