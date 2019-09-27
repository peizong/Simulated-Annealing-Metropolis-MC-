#!/bin/bash -l
#SBATCH -p regular
#SBATCH --qos=premium
#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH -L SCRATCH     #note: specify license need for the file systems your job needs, such as SCRATCH,project
#SBATCH -C haswell   #Use Haswell nodes
srun -n 1 ~/atat/src/wlce
