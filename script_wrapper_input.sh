#!/bin/bash
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=daisuke.shimaoka@monash.edu
#SBATCH --job-name=Wrapper_input_array
#SBATCH --time=40:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=33000
#SBATCH --array=1-4
module load matlab
matlab -nodisplay -nodesktop -nosplash < wrapper_input_array.m
