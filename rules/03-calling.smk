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
		ref = config['ref'],
		bam = "data/interim/add_rg/{sample}.rg.dedup.bam",
		bai = "data/interim/add_rg/{sample}.rg.dedup.bai"
	output:
		"data/interim/gvcf_files_bpres/{sample}.raw.snps.indels.g.vcf"
	params:
		regions = config['interval_list']
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
		ref = config['ref']
	output:
		expand("data/processed/scattered_intervals/{intervals}-scattered.intervals",
		intervals = INTERVALS)
	params:
		regions = config['interval_list'],
		intervals = len(INTERVALS),
		outdir = "data/processed/scattered_intervals"
	run:
		shell("gatk SplitIntervals \
		-R {input.ref} \
		-L {params.regions} \
		--scatter-count {params.intervals} \
		-O {params.outdir}")

################################################################################
# combine GVCFs with GenomicsDBImport
# https://software.broadinstitute.org/gatk/documentation/article?id=11813
# recommendation is to scatter over number of intervals ~ number of samples
# snakemake considerations:
# 	https://bitbucket.org/snakemake/snakemake/issues/895/combine-multiple-files-for-input-but
################################################################################

# files = lambda wildcards, input: " -V ".join(input),
# 		-V {params.files} \

rule combine_gvcfs:
	input:
		gvcfs = expand("data/interim/gvcf_files_bpres/{sample}.raw.snps.indels.g.vcf", sample = SAMPLES),
		region = "data/processed/scattered_intervals/{intervals}-scattered.intervals",
		map = "data/processed/sample_map"
	output:
		directory("data/interim/combined_database_bpres/{intervals}")
	params:
		tmp = "/scratch/sdturner/genomicsdbimport/{intervals}"
	run:
		shell("mkdir -p {params.tmp}")
		shell("gatk --java-options \"-Xmx60g -Xms60g\" \
		GenomicsDBImport \
		--genomicsdb-workspace-path {output} \
		--batch-size 50 \
		--reader-threads 6 \
		--sample-name-map {input.map} \
		--intervals {input.region} \
		--tmp-dir {params.tmp}")
		shell("rm -rf {params.tmp}")

################################################################################
# joint genotyping to produce VCF (raw SNPs & indels)
################################################################################

rule joint_geno:
	input:
		dir = "data/interim/combined_database_bpres/{intervals}",
		ref = config['ref']
	output:
		"data/raw/vcf_bpres/{intervals}.raw.snps.indels.vcf"
	params:
		db = "gendb://data/interim/combined_database_bpres/{intervals}",
		region = "data/processed/scattered_intervals/{intervals}-scattered.intervals"
	run:
		shell("gatk GenotypeGVCFs \
		-R {input.ref} \
		-V {params.db} \
		-L {params.region} \
		-new-qual \
		-G StandardAnnotation \
		-G AS_StandardAnnotation \
		--include-non-variant-sites \
		-O {output}")

### START HERE
rule vcf_concat:
	input:
		expand("data/raw/vcf_bpres/{intervals}.raw.snps.indels.vcf", intervals = INTERVALS)
	output:
		"data/interim/all_samples_unfiltered.vcf.gz"
	run:
		shell("bcftools concat {input} -Oz -o {output}")
		shell("tabix -p vcf {output}")
