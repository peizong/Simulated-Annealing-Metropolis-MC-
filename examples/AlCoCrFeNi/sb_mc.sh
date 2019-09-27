#!/bin/bash -l
#SBATCH --job-name="My Program"
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=general
##One node has two CPUs (2 tasks per node). One task will be placed on each CPU with 20 OpenMP threads
#module load intel/2019.1.053 intelmpi/2019.1.14

module load intel/2019.1.053 intelmpi/2019.1.144

export OMP_NUM_THREADS=40

#time  $HOME/atat/src/mcce-metropolis-test
time mpiexec.hydra -np 1 $HOME/atat/src/mcce-metropolis-new-swap
#time mpiexec.hydra -np 1 $HOME/atat/src/mcce-metropolis-test
