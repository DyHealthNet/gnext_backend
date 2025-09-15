#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=dyhealthnet_build_chr_bgz
#SBATCH --error=logging/%x_%j.err
#SBATCH --output=logging/%x.%j.log
#SBATCH --mem=64G
#SBATCH --array=1-22
#SBATCH --partition=slow-mc2

# Configure paths
OUT_DIR=/storage03/larend/pgwas_data/pgwas_out/GWAS_chr_bgz
mkdir -p "${OUT_DIR}"
CHR=$SLURM_ARRAY_TASK_ID

# Run the builder
srun python ../manage.py build_chr_bgz \
  --chrom "${CHR}" \
  --out-dir "${OUT_DIR}" \
  --row-block-size 64000 \
  --trait-workers 16
