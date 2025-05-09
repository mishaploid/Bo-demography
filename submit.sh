#!/bin/bash

date=$(date "+%Y_%m_%d")
# file="snakemake_logs/${date}_bodem.log"
echo $date
# echo $file

module load gsl/2.7.1
module load fastqc
module load trimmomatic/0.39
module load R
module load maven
module load openjdk/16.0.2
module load gatk
module load qualimap/2.3
module load vcftools
module load plink-ng/2.00a3.7
# module load conda/plink
# module load popvae
# module load singularity
module load apptainer
module load bedtools2/2.31.1

snakemake --jobs 200 --rerun-triggers mtime \
--use-conda \
--use-singularity \
--rerun-incomplete \
--latency-wait 60 \
--cluster-config submit.json \
--cluster "sbatch -v CONDA_AUTO_ACTIVATE_BASE=false -p {cluster.p} -o {cluster.o} --mem {cluster.mem} --time {cluster.time} --job-name {cluster.name}" # -p &>> $file
