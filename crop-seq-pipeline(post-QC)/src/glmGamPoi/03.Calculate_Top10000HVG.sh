#!/bin/bash

#SBATCH --job-name=Top10000Genes_Sublib2
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=200G
#SBATCH -t 2:00:00
#SBATCH --array=1,2,4 
#SBATCH -o /gstore/scratch/u/kharitoe/sublib1_2_4_log/Top10000Genes_Sublib_%A_%a.out        # <--- will create an output log file
#SBATCH -e /gstore/scratch/u/kharitoe/sublib1_2_4_log/Top10000Genes_Sublib_%A_%a.err        # <--- will create an error log file
#SBATCH --mail-type=END,FAIL       # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=kharitoe@gene.com 


#ml spaces/gpy
#ml gpy_gpu/Python310 


## Script to run
script=/gstore/project/elenakharitonova_melocars/downsampling_perturbseq/downsampling_github/glmGamPoi_Top10000_Genes/03.CalculateTop10000_genes_bysublib.py ## Calculate Top 10000 Genes
#script=/gstore/project/elenakharitonova_melocars/downsampling_perturbseq/downsampling_github/glmGamPoi_Top10000_Genes/03.CalculateTop10000_genes_concat.py ## Calculate Top 10000 Genes
#Note: Concat version of scripts requires no input

## Submit analysis
/apps/user/gpy/envs/gpy_gpu/Python310/bin/python3.10 -u ${script} -s sublib$SLURM_ARRAY_TASK_ID #Calculate Top 10000 HVG
