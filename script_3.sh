#!/bin/bash
#
#SBATCH -n 8
#SBATCH --partition=fast
#SBATCH --mem 20GB

module load python/3.9

python3 ext_script/script_3.py $1
