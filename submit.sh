#!/bin/bash

date=$(date "+%Y_%m_%d")
file="snakemake_logs/${date}_bodem.log"
echo $date
echo $file

module load trimmomatic/0.36
module load R
module load maven
module load java
module load GATK/4.0
module load bcftools
module load plink/1.90

snakemake --jobs 120 --use-conda \
--rerun-incomplete \
--latency-wait 60 \
--cluster-config submit.json \
--cluster "sbatch -p {cluster.p} -o {cluster.o} --cpus-per-task {cluster.cpus-per-task} --time {cluster.time}" -p &>> $file
