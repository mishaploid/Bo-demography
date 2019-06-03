## SNP calling pipeline for Brassica oleracea
## Follows GATK best practices https://software.broadinstitute.org/gatk/documentation/article.php?id=3893
## Sarah Turner-Hissong
## 4 March 2019

# requires a conda environment to run:
# conda env create --name bo-demography --file environment.yaml

# to run workflow use submit.sh script

# get sample names (option 1: pull unique identifier from fastq files)
# SAMPLES = glob_wildcards("data/external/fastq_raw/{sample}_pass_1.fastq.gz").sample

# get sample names (option 2: manual entries)
SAMPLES = ['SamC_' + str(x).rjust(3, '0') for x in range(1,120)]
SAMPLES2 = SAMPLES + ['B_rapa']

# print(SAMPLES)

# dictionary of SRA identifiers from downloaded project info
import csv

with open('data/external/Sra_oleracea.csv', mode = 'r') as infile:
	reader = csv.reader(infile)
	sample_dict = {rows[29]:rows[0] for rows in reader}

# add Brassica rapa reference (Chiifu)
sample_dict.update({'B_rapa':'SRR7881031'})

# list chromosome ids
chr = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']

# add paths for locally installed software
# ngsTools (https://github.com/mfumagalli/ngsTools)
ngsTools = "../software/ngsTools"
# FINESTRUCTURE software (v4)
# wget https://people.maths.bris.ac.uk/~madjl/finestructure/fs_4.0.1.zip
# unzip fs_4.0.1.zip
fs = "../bin/fs_4.0.1"

# a pseudo-rule that collects the target files (expected outputs)
rule all:
	input:
		# get_ref
		expand("data/external/ref/Boleracea_chromosomes.fasta"),
		# download SRA data, trim ends, and map to reference
		# expand("data/interim/mapped_reads/{sample}.bam", sample = SAMPLES2),
		# sort_bam - output is sorted BAM files for each sample
		expand("data/raw/sorted_reads/{sample}.sorted.bam", sample = SAMPLES2),
		# mark_dups - output is BAM files with duplicates marked
		# expand("data/interim/mark_dups/{sample}.dedup.bam", sample = SAMPLES),
		# add_rg - add read group information to BAM files
		# expand("data/interim/add_rg/{sample}.rg.dedup.bam", sample = SAMPLES),
		# coverage_depth
		"reports/bamqc/report.pdf",
		# hap_caller - output is GVCF file for each sample
		expand("data/interim/{sample}.raw.snps.indels.g.vcf", sample = SAMPLES)
		# combine_gvcfs - database for combining multiple gvcf files
		# expand("data/interim/combined_database/{chr}/vcfheader.vcf", chr = chr),
		# directory(expand("data/interim/combined_database/{chr}", chr = chr))
		# # joint_geno - outputs joint SNP calls for gvcf files
		# expand("data/raw/{chr}.raw.snps.indels.vcf", chr = chr),
		# # # get_snps
		# expand("data/processed/{chr}.filtered.snps.vcf", chr = chr),
		# # # bgzip_vcf
		# expand("data/processed/{chr}.filtered.snps.vcf.gz", chr = chr),
		# # # admix_input
		# "models/admixture/combined.pruned.bed",
		# # # admixture
		# expand("models/admixture/combined.pruned.{k}.Q", k = [1,2,3,4,5,7,8]),
		# # angsd_depth
		# expand("reports/ALL.{chr}.qc.depthGlobal", chr = chr)


include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/filtering.smk"
include: "rules/admixture.smk"
include: "rules/angsd.smk"
