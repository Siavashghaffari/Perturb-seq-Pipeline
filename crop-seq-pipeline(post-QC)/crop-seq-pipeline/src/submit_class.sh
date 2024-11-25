#!/bin/bash

#SBATCH -p himem
#SBATCH --mem=512G
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH -c 8


ml spaces/gpy
ml gpydev/gpuy310
#ml gpydev/gpy39
python ./Main_class.py