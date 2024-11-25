#!/bin/bash

#SBATCH -p himem
#SBATCH --mem=512G
#SBATCH --qos=long
#SBATCH -c 8


ml spaces/gpy
ml gpydev/gpuy310
#ml gpydev/gpy39
python ./Main_ED2.py