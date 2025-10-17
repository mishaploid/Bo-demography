#!/bin/bash

date=$(date "+%Y_%m_%d")
echo $date

# R and GATK are in module jail because they interfere with the default python3 path and that is RUDE
module load gsl/2.7.1
# module load fastqc
# module load trimmomatic/0.39
# module load R
# module load maven
module load openjdk/16.0.2
# module load gatk
# module load qualimap/2.3
# module load vcftools
# module load apptainer
# module load bedtools2/2.31.1
# module load plink-ng/2.00a3.7
# module load conda/plink
# module load popvae
# module load singularity

snakemake --profile slurm \
--workflow-profile ./workflow_profile \
--software-deployment-method conda apptainer \
--rerun-triggers mtime 
