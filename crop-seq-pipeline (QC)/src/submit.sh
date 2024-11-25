#!/bin/bash

#SBATCH -p himem
#SBATCH --mem=256G
#SBATCH -c 8


ml spaces/gpy
ml gpydev/gpuy310
#ml gpydev/gpy39
python ./Main.py