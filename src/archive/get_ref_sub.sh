#!/bin/bash -l
#SBATCH -D /home/turnersa/data/bo-demography/sequence_data/reference
#SBATCH -J get_ref
#SBATCH -o slurm-logs/get_ref%j.out
#SBATCH -e slurm-logs/get_ref%j.err

set -e
set -u

./get_ref.sh
