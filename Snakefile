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
SAMPLES = glob_wildcards("data/external/fastq_raw/{sample}_pass_1.fastq.gz").sample
SAMPLES2 = [s for s in SAMPLES if s not in {"SRR7881031"}]

print(SAMPLES)
print(SAMPLES2)

# chromosome ids
chr = ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"]

# paths to required software executables
trimmomatic = "../software/Trimmomatic-0.38/trimmomatic-0.38.jar"
gatk = "../software/gatk-4.0.12.0"

# expected outputs
rule all:
	input:
		# trimmomatic_pe - output is FASTQ files with adapters trimmed
		# expand("data/interim/trimmed_reads/{sample}.{pair}.fastq.gz", sample = SAMPLES, pair = [1,2]),
		# bwa_map - output is aligned BAM files for each sample
		# expand("data/interim/mapped_reads/{sample}.bam", sample = SAMPLES),
		# sort_bam - output is sorted BAM files for each sample
		expand("data/raw/sorted_reads/{sample}.bam", sample = SAMPLES),
		# mark_dups - output is BAM files with duplicates marked
		# expand("data/interim/dedup_reads/{sample}_dedup.bam", sample = SAMPLES),
		# add_rg - add read group information to BAM files
		# expand("data/interim/dedup_rg/{sample}_dedup_rg.bam", sample = SAMPLES),
		# hap_caller - output is GVCF file for each sample
		expand("data/interim/{sample}.raw.snps.indels.g.vcf", sample = SAMPLES2),
		# combine_gvcfs - database for combining multiple gvcf files
		directory("data/interim/combined_database"),
		# joint_geno - outputs joint SNP calls for gvcf files
		# "data/raw/raw_variants.vcf",
		# get_snps
		"data/raw/raw_snps.vcf",
		# filter_snps
		"data/processed/filtered_snps.vcf"

### Step 1: Data cleanup

# trim adapters
rule trimmomatic_pe:
    input:
        "data/external/fastq_raw/{sample}_pass_1.fastq.gz",
        "data/external/fastq_raw/{sample}_pass_2.fastq.gz"
    output:
        temp("data/interim/trimmed_reads/{sample}.1.fastq.gz"),
        temp("data/interim/trimmed_reads/{sample}.1.unpaired.fastq.gz"),
        temp("data/interim/trimmed_reads/{sample}.2.fastq.gz"),
        temp("data/interim/trimmed_reads/{sample}.2.unpaired.fastq.gz")
    log:
        "logs/trimmomatic/{sample}.log"
    run:
        shell("java -jar {trimmomatic} PE {input} {output} \
		LEADING:3 \
		TRAILING:3 \
		SLIDINGWINDOW:4:15 \
		MINLEN:36 2> {log}")

# LEADING:3 - remove leading low quality or N bases (quality < 3)
# TRAILING:3 - remove trailing low quality or N bases (quality < 3)
# SLIDINGWINDOW: scan read with 4-base sliding window, cutting when average quality
# drops below 15
# MINLEN:36 - drop reads below 36 bases long
# note: check sizes of paired vs. unpaired output files (expect paired to be larger)

# map to reference (TO1000 v2.1 from Parkin et al. 2014)
rule bwa_map:
    input:
       "data/external/ref/Brassica_oleracea.v2.1.dna.toplevel.fa",
       "data/interim/trimmed_reads/{sample}.1.fastq.gz",
       "data/interim/trimmed_reads/{sample}.2.fastq.gz"
    output:
       temp("data/interim/mapped_reads/{sample}.bam")
	log:
       "logs/bwa/{sample}.log"
	threads: 8
	run:
		shell("(bwa mem -t {threads} {input} | \
		samtools view -Sb > {output}) 2> {log}")

# sort BAM files with Picard
rule sort_bam:
	input:
		"data/interim/mapped_reads/{sample}.bam"
	output:
		"data/raw/sorted_reads/{sample}.bam"
	params:
		tmp = "data/interim/sorted_reads/tmp"
	run:
		shell("{gatk}/gatk SortSam \
		-I={input} \
		-O={output} \
		--SORT_ORDER=coordinate \
		--TMP_DIR={params.tmp}")

# mark duplicates with Picard
# no need to remove duplicates here - Haplotype Caller will ignore them
rule mark_dups:
	input:
		"data/raw/sorted_reads/{sample}.bam"
	output:
		bam = temp("data/interim/dedup_reads/{sample}_dedup.bam"),
		index = temp("data/interim/dedup_reads/{sample}_dedup.bai"),
		metrics = temp("data/interim/dedup_reads/{sample}_metrics.txt"),
	params:
		tmp = "data/interim/dedup_reads/tmp"
	run:
		shell("{gatk}/gatk MarkDuplicates \
		-I={input} \
		-O={output.bam} \
		--METRICS_FILE={output.metrics} \
		--CREATE_INDEX=true \
		-MAX_FILE_HANDLES=1000 \
		--ASSUME_SORT_ORDER=coordinate \
		--TMP_DIR={params.tmp}")
# metrics file - records duplication metrics (required)
# create_index - true; indexes resulting BAM file
# -MAX_FILE_HANDLES - reduced this to meet system requirements for max number of open files (check with `ulimit -n`)

# add read groups with Picard
# info is just dummy variables
rule add_rg:
	input:
		"data/interim/dedup_reads/{sample}_dedup.bam"
	output:
		bam = temp("data/interim/dedup_rg/{sample}_dedup_rg.bam"),
		index = temp("data/interim/dedup_rg/{sample}_dedup_rg.bai")
	params:
		tmp = "data/raw/bam_files/tmp"
	run:
		shell("{gatk}/gatk AddOrReplaceReadGroups \
		-I={input} \
		-O={output.bam} \
		--CREATE_INDEX=true \
		-RGID=4 \
		-RGLB=lib1 \
		-RGPL=illumina \
		-RGPU=unit1 \
		-RGSM=20 \
		--TMP_DIR {params.tmp}")

# Step 2: variant discovery!

# haplotype caller
# https://software.broadinstitute.org/gatk/documentation/tooldocs/3.6-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php

rule hap_caller:
	input:
		ref = "data/external/ref/Brassica_oleracea.v2.1.dna.toplevel.fa",
		bam = "data/interim/dedup_rg/{sample}_dedup_rg.bam"
	output:
		"data/interim/{sample}.raw.snps.indels.g.vcf"
	run:
		shell("{gatk}/gatk HaplotypeCaller \
		-I {input.bam} \
		-O {output} \
		-R {input.ref} \
		-G StandardAnnotation \
		-G AS_StandardAnnotation \
		-G StandardHCAnnotation \
		--emit-ref-confidence GVCF")

# combine GVCFs
# https://software.broadinstitute.org/gatk/documentation/article?id=11813
# snakemake considerations - https://bitbucket.org/snakemake/snakemake/issues/895/combine-multiple-files-for-input-but
# expand("data/interim/{sample}.raw.snps.indels.g.vcf", sample = SAMPLES.remove("SRR7881031"))
rule combine_gvcfs:
	input:
		expand("data/interim/{sample}.raw.snps.indels.g.vcf", sample = SAMPLES2)
	output:
		directory("data/interim/combined_database")
	params:
		files = lambda wildcards, input: " -V ".join(input)
	run:
		shell("{gatk}/gatk GenomicsDBImport \
		-V {params.files} \
		-G StandardAnnotation \
		-G AS_StandardAnnotation \
		--genomicsdb-workspace-path {output} \
		--intervals C1,C2,C3,C4,C5,C6,C7,C8,C9")

# joint genotyping - raw SNP and indel VCF
rule joint_geno:
	input:
		ref = "data/external/ref/Brassica_oleracea.v2.1.dna.toplevel.fa"
	output:
		"data/raw/raw_variants.vcf"
	params:
		db = "gendb://data/interim/combined_database"
	run:
		shell("{gatk}/gatk GenotypeGVCFs \
		-R {input.ref} \
		-V {params.db} \
		-newQual \
		-G StandardAnnotation \
		-G AS_StandardAnnotation \
		-O {output}")

# hard filtering
# https://software.broadinstitute.org/gatk/documentation/article.php?id=2806
# Note: Variant Quality Score Recalibration is not an option (need truth/training sets)

# extract SNPs
rule get_snps:
	input:
		ref = "data/external/ref/Brassica_oleracea.v2.1.dna.toplevel.fa",
		vcf = "data/raw/raw_variants.vcf"
	output:
		"data/raw/raw_snps.vcf"
	run:
		shell("{gatk}/gatk SelectVariants \
		-R {input.ref} \
		-V {input.vcf} \
		-select-type SNP \
		-O {output}")

# apply filters
# using base GATK recommendations
# https://software.broadinstitute.org/gatk/documentation/article?id=23216#2

rule filter_snps:
	input:
		ref = "data/external/ref/Brassica_oleracea.v2.1.dna.toplevel.fa",
		vcf = "data/raw/raw_snps.vcf"
	output:
		"data/processed/filtered_snps.vcf"
	run:
		shell("{gatk}/gatk VariantFiltration \
		-filter \"QD < 2.0\" --filter-name \"QD2\" \
		-filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \
		-filter \"SOR > 3.0\" --filter-name \"SOR3\" \
		-filter \"FS > 60.0\" --filter-name \"FS60\" \
		-filter \"MQ < 40.0\" --filter-name \"MQ40\" \
		-filter \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" \
		-filter \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8\" \
		-O {output}")

# evaluate filtering parameters
# https://software.broadinstitute.org/gatk/documentation/article.php?id=11069
