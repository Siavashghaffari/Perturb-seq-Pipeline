#!/bin/bash

# SBATCH --job-name=Bin_Dev
# SBATCH -n 10
# SBATCH -N 1
# SBATCH --mem=700G
# SBATCH --qos=short
# SBATCH --partition himem
# SBATCH --mail-type=FAIL       # Type of email notification- BEGIN,END,FAIL,ALL
# SBATCH --mail-user=ghaffars@gene.com 


ml R/dev



# /gstore/home/ghaffars/Projects/bdev/cropquestwizard/Rcane/wrappers/CommandLineMPE.R --numGenesToKeep 3000 --nworkers 10 --sceName #/gstore/scratch/u/ghaffars/Dataset/sublib1/raw_qc.h5ad

/gstore/home/ghaffars/Projects/bdev/cropquestwizard/Rcane/wrappers/CommandLineMPE.R --numGenesToKeep 3000 --nworkers 10 --sceName "DS000017114"


