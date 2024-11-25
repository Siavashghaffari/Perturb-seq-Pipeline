#!/bin/bash

#SBATCH --job-name=Shrink_Betas
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --mail-type=FAIL      
#SBATCH --mail-user=ghaffars@gene.com 

module load R/dev


################################################################################################################
########################################### Specify Inputs #####################################################
################################################################################################################


## Script to run
script=/gstore/home/ghaffars/Cumulus/crc_dld1_sublib1_bdev/crop-seq-pipeline/src/glmGamPoi_NE/07.Shrink_Beta_Estimates.R



################################################################################################################
###################################### Run Scripts #############################################################
################################################################################################################


## Run analysis
Rscript ${script} -s sublib1 -niter 32 ## Run for sublibrary 1

