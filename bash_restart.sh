#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gpus-per-node=4
#SBATCH --ntasks-per-node=4            # ranks per node, so 4 means each gpu gets a rank
#SBATCH --partition=compute_full_node
#SBATCH --time=04:00:00                 # hours, minutes, seconds

#SBATCH -o outjob_test.o%j       # Name of stdout output file
#SBATCH -e outjob_test.e%j       # Name of stderr error file
#SBATCH -J test

module load gcc cmake cuda openmpi
srun --mpi=pmix ./athena \
    -r ./rst/rst1_gamma1_67_beta1_S10000_gpa13_3_nx3200.00005.rst \
    -i rst1_gamma1_67_beta1_S10000_gpa13_3_nx3200.athinput