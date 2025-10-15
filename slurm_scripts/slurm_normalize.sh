#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name=dyhealthnet_normalize
#SBATCH --error=logging/%x_%j.err
#SBATCH --output=logging/%x.%j.log
#SBATCH --mem=64G
#SBATCH --partition=slow-mc2

python ../manage.py normalize