#!/bin/bash
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=20G
#SBATCH -t 02:30:00

#SBATCH --array=0-12

# Array job to compute all ingredients for the
# Fisher forecast.

date

module purge
module load gsl/2.6
module load fftw/gcc/3.3.4

cd $HOME/hmpdf/forecast

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:..
export OMP_NUM_THREADS=4

./main_forecast $SLURM_ARRAY_TASK_ID

date
