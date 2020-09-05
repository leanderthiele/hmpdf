#!/bin/bash
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=20G
#SBATCH -t 02:30:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=lthiele@princeton.edu

#SBATCH --array=0-13

date

module purge
module load gsl/2.6
module load fftw/gcc/3.3.4

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/lthiele/hmpdf
export OMP_NUM_THREADS=4

cd /home/lthiele/hmpdf/forecast

./main_forecast $SLURM_ARRAY_TASK_ID

date
