#!/bin/bash

date=$(date "+%Y_%m_%d")
# file="snakemake_logs/${date}_bodem.log"
echo $date
# echo $file

module load fastqc
module load trimmomatic/0.39
module load R
module load maven
module load java
module load gatk
module load bcftools
module load qualimap/2.1.1
module load vcftools
module load plink/1.90
module load plink2
module load angsd
# module load popvae
module load singularity

snakemake --jobs 200 --use-conda \
--use-singularity \
--rerun-incomplete \
--latency-wait 60 \
--cluster-config submit.json \
--cluster "sbatch -p {cluster.p} -o {cluster.o} --mem {cluster.mem} --time {cluster.time} --job-name {cluster.name}" # -p &>> $file
