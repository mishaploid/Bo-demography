## SNP calling pipeline for Brassica oleracea
## Follows GATK best practices https://software.broadinstitute.org/gatk/documentation/article.php?id=3893
## Sarah Turner-Hissong
## 4 March 2019

# requires a conda environment to run:
# conda env create --name bo-demography --file environment.yaml

# to run full script:
# snakemake --use-conda --jobs N --rerun-incomplete --cluster-config submit.json --cluster "sbatch"
# a pseudo-rule that collects the target files

# get sample names
# SAMPLES = glob_wildcards("data/external/fastq_raw/{sample}_pass_1.fastq.gz").sample

# SAMPLES = ['SamC_' + str(x).rjust(3, '0') for x in range(1,120)]
# SAMPLES.append('SRR7881031')
# SAMPLES = ['SamC_010', 'SamC_011', 'SamC_013', 'SamC_024', 'SamC_028']
# SAMPLES2 = [s for s in SAMPLES if s not in {"SRR7881031"}]


SAMPLES = ['SamC_' + str(x).rjust(3, '0') for x in range(1,120)]

print(SAMPLES)

# chromosome ids
# chr = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
chr = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
# paths to required software executables
# trimmomatic = "../software/Trimmomatic-0.38/trimmomatic-0.38.jar"
# gatk = "../software/gatk-4.0.12.0"

# expected outputs
rule all:
	input:
		# trimmomatic_pe - output is FASTQ files with adapters trimmed
		# expand("data/interim/trimmed_reads/{sample}.{pair}.fastq.gz", sample = SAMPLES, pair = [1,2]),
		# bwa_map - output is aligned BAM files for each sample
		# expand("data/interim/mapped_reads/{sample}.bam", sample = SAMPLES),
		# sort_bam - output is sorted BAM files for each sample
		# expand("data/raw/sorted_reads/{sample}.sorted.bam", sample = SAMPLES),
		# mark_dups - output is BAM files with duplicates marked
		# expand("data/interim/dedup/{sample}.dedup.bam", sample = SAMPLES),
		# add_rg - add read group information to BAM files
		# expand("data/interim/add_rg/{sample}.rg.dedup.bam", sample = SAMPLES),
		# hap_caller - output is GVCF file for each sample
		# expand("data/interim/{sample}.raw.snps.indels.g.vcf", sample = SAMPLES),
		# rename samples
		expand("data/interim/{sample}.renamed.raw.snps.indels.g.vcf", sample = SAMPLES),
		# combine_gvcfs - database for combining multiple gvcf files
		directory(expand("data/interim/combined_database/{chr}", chr = chr)),
		# joint_geno - outputs joint SNP calls for gvcf files
		expand("data/raw/{chr}.raw.snps.indels.vcf", chr = chr),
		# get_snps
		expand("data/processed/{chr}.filtered.snps.vcf", chr = chr)
		# # filter_snps
		# "data/processed/filtered_snps.vcf"

include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/filtering.smk"
