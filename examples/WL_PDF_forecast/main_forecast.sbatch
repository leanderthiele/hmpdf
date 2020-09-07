#!/bin/bash

# Array job to compute all ingredients for the Fisher forecast.
# The requested resources are relatively conservative.

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=10G
#SBATCH -t 01:00:00
#SBATCH --array=0-12


date

module purge
module load gsl/2.6
module load fftw/gcc/3.3.4

# make libhmpdf.so visible to the linker
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../..

# for CLASS multithreading
export OMP_NUM_THREADS=4

# create the output directory
mkdir -p results

# run the executable for this index
./main_forecast $SLURM_ARRAY_TASK_ID

date
