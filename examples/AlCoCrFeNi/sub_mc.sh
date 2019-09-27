#!/bin/csh
## Change into the current working directory
#$ -cwd
##
## The name for the job. It will be displayed this way on qstat
#$ -N aimd
##
## Number of cores to request
#$ -pe mpi 1
##
#$ -r n
##
## Queue Name
#$ -q alb

##Remember to update the queue name!

##Load Modules
#module load intel/2014.1.046 intelmpi/4.1.3.048 vasp/5.4.1_gamma
#module load  intelmpi/2019.0.117
#module load  intel/2019.0.045
#module load  vasp/5.4.4 
##Run the parallel job
#mpirun vasp > output

time $HOME/atat/src/mcce-metropolis
