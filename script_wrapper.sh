#!/bin/bash
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=daisuke.shimaoka@monash.edu
#SBATCH --job-name=Wrapper
#SBATCH --time=20:00:00
#SBATCH --ntasks=1-81
#SBATCH --mem-per-cpu=33000
module load matlab
matlab -nodisplay -nodesktop -nosplash < wrapper_array.m
