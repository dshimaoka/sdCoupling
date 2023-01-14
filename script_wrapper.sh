#!/bin/bash
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=daisuke.shimaoka@monash.edu
#SBATCH --job-name=Wrapper_array
#SBATCH --time=40:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=33000
#SBATCH --array=1-10
module load matlab
matlab -nodisplay -nodesktop -nosplash < wrapper_array.m
