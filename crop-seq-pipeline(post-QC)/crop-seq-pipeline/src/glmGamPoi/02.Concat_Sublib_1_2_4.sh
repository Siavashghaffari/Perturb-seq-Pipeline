#!/bin/bash

#SBATCH --job-name=Concate_Sublibs
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p cpu
#SBATCH --mem=500G
#SBATCH --qos=short
#SBATCH -o /gstore/scratch/u/kharitoe/sublib1_2_4_log/Concat_Sublibs_%A_%a.out        # <--- will create an output log file
#SBATCH -e /gstore/scratch/u/kharitoe/sublib1_2_4_log/Concat_Sublibs_%A_%a.err        # <--- will create an error log file
#SBATCH --mail-type=END,FAIL       # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=kharitoe@gene.com 


#ml spaces/gpy
#ml gpy_cpu/Python310 


## Script to run
script1=/gstore/project/elenakharitonova_melocars/downsampling_perturbseq/downsampling_github/glmGamPoi_Top10000_Genes/02.Concat_Sublib_1_2_4.py ## Calculate Top 2000 Genes

# Exit immediately if a command exits with a non-zero status.
set -e

/apps/user/gpy/envs/gpy_cpu/Python310/bin/python3.10 -u ${script1}  # Calculate Top 2000 Genes
   