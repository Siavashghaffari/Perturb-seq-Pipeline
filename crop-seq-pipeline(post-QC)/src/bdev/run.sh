#!/bin/bash

#SBATCH -p himem
#SBATCH --mem=512G
#SBATCH -c 8

ml spaces/gpy
ml gpy_gpu/Python310

python ./split.py