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
SAMPLES = ['SamC_010']
SAMPLES2 = [s for s in SAMPLES if s not in {"SRR7881031"}]

print(SAMPLES)
print(SAMPLES2)

# chromosome ids
chr = ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"]

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
		expand("data/raw/sorted_reads/{sample}.sorted.bam", sample = SAMPLES),
		# mark_dups - output is BAM files with duplicates marked
		expand(temp("data/raw/{sample}.dedup.bam"), sample = SAMPLES),
		# add_rg - add read group information to BAM files
		expand(temp("data/raw/{sample}.rg.dedup.bam"), sample = SAMPLES),
		# hap_caller - output is GVCF file for each sample
		expand("data/sdturner/{sample}.raw.snps.indels.g.vcf", sample = SAMPLES)
		# combine_gvcfs - database for combining multiple gvcf files
		# directory("data/interim/combined_database"),
		# joint_geno - outputs joint SNP calls for gvcf files
		# "data/raw/raw_variants.vcf"
		# get_snps
		# "data/raw/raw_snps.vcf",
		# filter_snps
		# "data/processed/filtered_snps.vcf"

include: "rules/mapping.smk"
include: "rules/calling.smk"




# hard filtering
# https://software.broadinstitute.org/gatk/documentation/article.php?id=2806
# Note: Variant Quality Score Recalibration is not an option (need truth/training sets)

# extract SNPs
# rule get_snps:
# 	input:
# 		ref = "data/external/ref/Brassica_oleracea.v2.1.dna.toplevel.fa",
# 		vcf = "data/raw/raw_variants.vcf"
# 	output:
# 		temp("data/raw/raw_snps.vcf")
# 	run:
# 		shell("{gatk}/gatk SelectVariants \
# 		-R {input.ref} \
# 		-V {input.vcf} \
# 		-select-type SNP \
# 		-O {output}")
#
# # apply filters
# # using base GATK recommendations
# # https://software.broadinstitute.org/gatk/documentation/article?id=23216#2
#
# rule filter_snps:
# 	input:
# 		ref = "data/external/ref/Brassica_oleracea.v2.1.dna.toplevel.fa",
# 		vcf = "data/raw/raw_snps.vcf"
# 	output:
# 		"data/processed/filtered_snps.vcf"
# 	run:
# 		shell("{gatk}/gatk VariantFiltration \
# 		-filter \"QD < 2.0\" --filter-name \"QD2\" \
# 		-filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \
# 		-filter \"SOR > 3.0\" --filter-name \"SOR3\" \
# 		-filter \"FS > 60.0\" --filter-name \"FS60\" \
# 		-filter \"MQ < 40.0\" --filter-name \"MQ40\" \
# 		-filter \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" \
# 		-filter \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8\" \
# 		-O {output}")
#
# # evaluate filtering parameters
# # https://software.broadinstitute.org/gatk/documentation/article.php?id=11069
