#!/bin/bash

#SBATCH --job-name=Run_bdev
#SBATCH -n 10
#SBATCH -N 1
#SBATCH --array=1-5000 ## Ideally set to be 1-5200 but interns can't submit job arrays of 5000, so I submit them separately
#SBATCH --mem=20G
#SBATCH --qos=short
#SBATCH -t 1-     
#SBATCH --mail-type=FAIL      
#SBATCH --mail-user=ghaffars@gene.com 


ml R/dev

################################################################################################################
########################################### Specify Inputs #####################################################
################################################################################################################





## Input file split by arget
input_directory=/gstore/scratch/u/ghaffars/glmGamPoi/sublib1_bdev/data/

## Gene File - Info of Gene Names - Aka Rows/Features 
gene_file=/gstore/scratch/u/ghaffars/glmGamPoi/sublib1_bdev/data/remove_pos_cont_gene_data.csv

## Output Directory where to store results
output_dir=/gstore/scratch/u/ghaffars/glmGamPoi/bdev/sublib1/



## List of targets (within each library)
target_list=$(tail ${input_directory}coldata_for_glmgampois.csv -n +2 | awk -F "," '{print $7}' | sort | uniq)
target_arr=($target_list)

## Pring which target you are running
echo ${SLURM_ARRAY_TASK_ID}
target=${target_arr[${SLURM_ARRAY_TASK_ID}-1]}
echo $target


################################################################################################################
###############                         Submit Script                       ##################################
################################################################################################################

/gstore/home/ghaffars/Projects/bdev/cropquestwizard/Rcane/wrappers/CommandLineMPE.R --numGenesToKeep 3000 --nworkers 10 --sceName ${input_directory}by_target/${target}.h5ad --binomialDevianceOutput ${output_dir}bdev_${target}.csv


