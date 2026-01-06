#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=4
#SBATCH --ntasks-per-node=4            # ranks per node, so 4 means each gpu gets a rank
#SBATCH --partition=compute_full_node
#SBATCH --time=02:00:00                 # hours, minutes, seconds

#SBATCH -o outjob_test.o%j       # Name of stdout output file
#SBATCH -e outjob_test.e%j       # Name of stderr error file
#SBATCH -J test

module load gcc cmake cuda openmpi
srun --mpi=pmix ./athena -i gamma1_50_beta0_001_S10000_gpa13_3_nx1600.athinput
srun --mpi=pmix ./athena -i gamma1_50_beta1000_S10000_gpa13_3_nx1600.athinput
srun --mpi=pmix ./athena -i gamma1_50_beta1_S10000_gpa13_3_nx1600.athinput