################################################################################
## Rules for variant discovery
## HaplotypeCaller -> GenomicsDBImport -> GenotypeGVCFs
################################################################################

################################################################################
# Run HaplotypeCaller for each sample
# calls SNPs and indels via local re-assembly of haplotypes
# able to call difficult regions
# -L = params file that specifies regions to call
# -G = annotations to include
#		StandardAnnotation
#		AS_StandardAnnotation (allele specific)
# --emit-ref-confidence = mode for emitting reference confidence scores (BP resolution chosen here)
# 		see details: https://software.broadinstitute.org/gatk/documentation/article.php?id=4017
################################################################################

rule hap_caller:
	input:
		ref = "data/external/ref/Boleracea_chromosomes.fasta",
		bam = "data/interim/add_rg/{sample}.rg.dedup.bam",
		bai = "data/interim/add_rg/{sample}.rg.dedup.bai"
	output:
		"data/interim/gvcf_files_bpres/{sample}.raw.snps.indels.g.vcf"
	params:
		regions = "data/raw/b_oleracea.interval_list"
	run:
		shell("gatk HaplotypeCaller \
		-I {input.bam} \
		-O {output} \
		-R {input.ref} \
		-L {params.regions} \
		-G StandardAnnotation \
		-G AS_StandardAnnotation \
		--emit-ref-confidence BP_RESOLUTION")

################################################################################
# scatter reference into intervals using SplitIntervals
# https://gatk.broadinstitute.org/hc/en-us/articles/360036348372-SplitIntervals
# ``--scatter-count n` splits reference into n intervals
################################################################################

rule split_intervals:
	input:
		ref = "data/external/ref/Boleracea_chromosomes.fasta"
	output:
		expand("data/processed/scattered_intervals/{count}-scattered.intervals",
		count = INTERVALS)
	params:
		regions = "data/raw/b_oleracea.interval_list",
		count = len(INTERVALS),
		outdir = "data/processed/scattered_intervals"
	run:
		shell("gatk SplitIntervals \
		-R {input.ref} \
		-L {params.regions} \
		--scatter-count {params.count} \
		-O {params.outdir}")

################################################################################
# combine GVCFs with GenomicsDBImport
# https://software.broadinstitute.org/gatk/documentation/article?id=11813
# recommendation is to scatter over number of intervals ~ number of samples
# snakemake considerations:
# 	https://bitbucket.org/snakemake/snakemake/issues/895/combine-multiple-files-for-input-but
################################################################################

rule combine_gvcfs:
	input:
		expand("data/interim/gvcf_files_bpres/{sample}.raw.snps.indels.g.vcf", sample = SAMPLES)
	output:
		directory("data/interim/combined_database_bpres/{count}")
	params:
		files = lambda wildcards, input: " -V ".join(input),
		region = "data/processed/scattered_intervals/{count}-scattered.intervals",
		tmp = "/scratch/sdturner/genomicsdbimport/{count}",
		map = "data/processed/sample_map"
	run:
		shell("mkdir -p {params.tmp}")
		shell("gatk GenomicsDBImport \
		--genomicsdb-workspace-path {output} \
		--batch-size 50 \
		--reader-threads 2 \
		--sample-name-map {params.map} \
		--intervals {params.region} \
		--tmp-dir {params.tmp}")
		shell("rm -rf {params.tmp}")

# 		-V {params.files} \

################################################################################
# joint genotyping to produce VCF (raw SNPs & indels)

rule joint_geno:
	input:
		dir = directory("data/interim/combined_database_bpres/{count}"),
		ref = "data/external/ref/Boleracea_chromosomes.fasta"
	output:
		"data/raw/vcf_bpres/{count}.raw.snps.indels.vcf"
	params:
		db = "gendb://data/interim/combined_database_bpres/{count}"
	run:
		shell("gatk GenotypeGVCFs \
		-R {input.ref} \
		-V {params.db} \
		-new-qual \
		-G StandardAnnotation \
		-G AS_StandardAnnotation \
		--include-non-variant-sites \
		-O {output}")
