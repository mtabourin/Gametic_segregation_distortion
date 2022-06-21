#!/bin/bash
#
#SBATCH -n 8
#SBATCH --partition=fast
#SBATCH --mem 20GB

module load r/4.1.1
module load python/3.9

python3 script_4.py $1
