#! /usr/bin/env bash
#PBS -j oe
#PBS -N FeCo
#PBS -l walltime=2:00:00
#PBS -l nodes=1
#PBS -A MAT020

#-------------------------------------------

#--Change directory to the working directory.

cd $PBS_O_WORKDIR

#module swap PrgEnv-intel PrgEnv-gnu
module swap PrgEnv-pgi PrgEnv-gnu

export OMP_NUM_THREADS=8

#--Run the executable.
date
time $HOME/atat/atat/src/wlce-sro
date

