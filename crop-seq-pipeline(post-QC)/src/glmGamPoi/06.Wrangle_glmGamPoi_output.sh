#!/bin/bash

#SBATCH --job-name=Wrangle_glmGamPoi
#SBATCH -n 10
#SBATCH -N 1
#SBATCH --mem=50G
#SBATCH --mail-type=END,FAIL       # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ghaffars@gene.com 


module load R/dev

################################################################################################################
########################################### Specify Inputs #####################################################
################################################################################################################

## Script to run
script=/gstore/home/ghaffars/Cumulus/crc_dld1_sublib1_bdev/crop-seq-pipeline/src/glmGamPoi/06.Wrangle_glmGamPoi_output.R

## Input file and output directory
targ_list_dir_sublib1=/gstore/scratch/u/ghaffars/glmGamPoi/sublib1_bdev/data/



################################################################################################################
###################################### Run Scripts #############################################################
################################################################################################################


## Run analysis
Rscript ${script} -tl ${targ_list_dir_sublib1} -s sublib1 ## Run for sublibrary 3

