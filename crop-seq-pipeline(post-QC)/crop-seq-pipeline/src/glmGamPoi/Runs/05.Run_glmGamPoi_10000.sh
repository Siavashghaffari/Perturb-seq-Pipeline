#!/bin/bash

#SBATCH --job-name=Run_glmGamPoi
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --array=1-2000 ## Ideally set to be 1-5200 but interns can't submit job arrays of 5000, so I submit them separately
#SBATCH --mem=20G
#SBATCH -t 1-     
#SBATCH --mail-type=FAIL      
#SBATCH --mail-user=ghaffars@gene.com 


module load R/dev

################################################################################################################
########################################### Specify Inputs #####################################################
################################################################################################################


## Script to run
script=/gstore/home/ghaffars/Cumulus/crc_dld1_sublib1_bdev/crop-seq-pipeline/src/glmGamPoi/05.Run_glmGamPoi_10000.R



## Input file split by arget
input_directory=/gstore/scratch/u/ghaffars/glmGamPoi/sublib1_bdev/data/

## Gene File - Info of Gene Names - Aka Rows/Features 
gene_file=/gstore/scratch/u/ghaffars/glmGamPoi/sublib1_bdev/data/remove_pos_cont_gene_data.csv

## Output Directory where to store results

output_dir=/gstore/scratch/u/ghaffars/glmGamPoi/sublib1_bdev/results/

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

## Check that output file exists, then begin creation. This is useful for resubmitting those jobs that failed. 
if [ ! -e ${output_dir}unweightedglmGamPoi_gembatch_results/glm_unweighted_target_${target}_top10000_sublib124union_genes.rds ]; then
    echo file not exist begin creation
    
    ## Submit job
    Rscript ${script} -i ${input_directory} -g ${gene_file} -o ${output_dir} -t ${target}

    else
    echo file exists. Do not run script to avoid rereation
    fi

