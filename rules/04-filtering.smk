# TODO: try out piping in snakemake to avoid writing files to disk
# https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#piped-output
# hard filtering
# https://software.broadinstitute.org/gatk/documentation/article.php?id=2806
# Note: Variant Quality Score Recalibration is not an option (need truth/training sets)

# remove indels
# select variant and invariant sites
# restrict alleles to be biallelic
rule select_variant:
	input:
		ref = config['ref'],
		vcf = "data/interim/all_samples_unfiltered.vcf.gz"
	output:
		"data/raw/vcf_bpres/all_samps.raw.snps.vcf"
	run:
		shell("gatk SelectVariants \
		-R {input.ref} \
		-V {input.vcf} \
		-select-type SNP \
		--restrict-alleles-to BIALLELIC \
		-O {output}")

rule select_invariant:
	input:
		ref = config['ref'],
		vcf = "data/interim/all_samples_unfiltered.vcf.gz"
	output:
		"data/raw/vcf_bpres/all_samps.raw.invariant.vcf"
	run:
		shell("gatk SelectVariants \
		-R {input.ref} \
		-V {input.vcf} \
		-select-type NO_VARIATION \
		-O {output}")

# filtering diagnostics
# https://evodify.com/gatk-in-non-model-organism/
# NOTE: some annotations cannot be calculated for invariant sites, because
# they require a mix of ref and alt reads
# Also, there is no GQ reported for invariant sites
# https://gatkforums.broadinstitute.org/gatk/discussion/4711/genotype-qualities-qual-scores-for-all-sites

rule vcf_diagnostics:
	input:
		ref = config['ref'],
		variant_vcf = "data/raw/vcf_bpres/all_samps.raw.snps.vcf",
		invariant_vcf = "data/raw/vcf_bpres/all_samps.raw.invariant.vcf"
	output:
		variant_stats = "reports/filtering_bpres/all_samps.raw.snps.table",
		invariant_stats = "reports/filtering_bpres/all_samps.raw.invariant.table",
		stats_plot = "reports/filtering_bpres/unfiltered_diagnostics.png"
	run:
		shell("gatk VariantsToTable \
		-R {input.ref} \
		-V {input.variant_vcf} \
		-F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F AD \
		-O {output.variant_stats}")
		shell("gatk VariantsToTable \
		-R {input.ref} \
		-V {input.invariant_vcf} \
		-F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F AD \
		-O {output.invariant_stats}")
		shell("Rscript src/filtering_diagnostics.R")


# Annotate variants with filtering criteria
# using base GATK recommendations

# Resources:
# https://software.broadinstitute.org/gatk/documentation/article?id=23216#2
# https://software.broadinstitute.org/gatk/documentation/article.php?id=3225
# https://informatics.fas.harvard.edu/whole-genome-resquencing-for-population-genomics-fastq-to-vcf.html#filtering
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531012--How-to-Filter-on-genotype-using-VariantFiltration

# Filters:
# QD = quality by depth
#	variant quality score divided by depth of alternate allele
# QUAL = SNP quality
# 	phred-scaled quality score; measure in confidence that there is a variation at a site
# SOR = StrandOddsRatio
#	high values suggest strand bias
# FS = Fisher strand
#	phred-scaled p-values for strand bias; higher values more likely to be false +
# MQ = RMS Mapping quality
#	root mean square of mapping quality of reads across samples
# MQRankSum = Mapping quality rank sum test
#	U-based z-approx from Mann-Whitney Rank Sum Test for MQ
# ReadPosRankSum
#	U-based z-approx from Mann-Whitney Rank Sum Test for distance of alt allele from end of reads
# 	bases near end of read more likely to have errors
# 	if all reads with the alt allele are near the end may be errors

rule filter_flags_snps:
	input:
		raw_snp_vcf = "data/raw/vcf_bpres/all_samps.raw.snps.vcf"
	output:
		snps_flagged = temp("data/interim/filtering/{chr}_allsamps.filtered.qual.dp6_200.snpsOnly.vcf"),
		snps_idx = temp("data/interim/filtering/{chr}_allsamps.filtered.qual.dp6_200.snpsOnly.vcf.idx")
	params:
		chr = "{chr}"
	run:
		shell("gatk VariantFiltration \
		-V {input.raw_snp_vcf} \
		--intervals {params.chr} \
		-filter \"QD < 2.0\" --filter-name \"QD2\" \
		-filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \
		-filter \"SOR > 3.0\" --filter-name \"SOR3\" \
		-filter \"FS > 60.0\" --filter-name \"FS60\" \
		-filter \"MQ < 40.0\" --filter-name \"MQ40\" \
		-filter \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" \
		-filter \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8\" \
		-G-filter \"DP < 6 || DP > 200\" \
		-G-filter-name \"DP_6-200\" \
		-O {output.snps_flagged}")

rule filtered_to_nocall_snps:
	input:
		snps_flagged = "data/interim/filtering/{chr}_allsamps.filtered.qual.dp6_200.snpsOnly.vcf",
		snps_idx = "data/interim/filtering/{chr}_allsamps.filtered.qual.dp6_200.snpsOnly.vcf.idx"
	output:
		snps_nocall = temp("data/interim/filtering/{chr}_allsamps.filtered.qual.dp6_200.snpsOnly.nocall.vcf"),
		snps_nocall_idx = temp("data/interim/filtering/{chr}_allsamps.filtered.qual.dp6_200.snpsOnly.nocall.vcf.idx")
	run:
		shell("gatk VariantFiltration \
		-V {input.snps_flagged} \
		--set-filtered-genotype-to-no-call true \
		-O {output.snps_nocall}")

rule filter_snps:
	input:
		snps_nocall = "data/interim/filtering/{chr}_allsamps.filtered.qual.dp6_200.snpsOnly.nocall.vcf",
		snps_nocall_idx = "data/interim/filtering/{chr}_allsamps.filtered.qual.dp6_200.snpsOnly.nocall.vcf.idx",
		ref = config['ref']
	output:
		filtered_snps = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp6_200.maxnocall10.biallelic.snps.vcf",
		filtered_snps_idx = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp6_200.maxnocall10.biallelic.snps.vcf.idx"
	run:
		shell("gatk SelectVariants \
		-R {input.ref} \
		-V {input.snps_nocall} \
		-select-type SNP \
		--restrict-alleles-to BIALLELIC \
		--max-nocall-fraction 0.1 \
		--exclude-filtered true \
		-O {output.filtered_snps}")

rule filter_flags_invariant:
	input:
		raw_invariant_vcf = "data/raw/vcf_bpres/all_samps.raw.invariant.vcf"
	output:
		invariant_flagged = temp("data/interim/filtering/{chr}_allsamps.filtered.qual.dp6_200.invariantOnly.vcf"),
		invariant_idx = temp("data/interim/filtering/{chr}_allsamps.filtered.qual.dp6_200.invariantOnly.vcf.idx")
	params:
		chr = "{chr}"
	run:
		shell("gatk VariantFiltration \
		-V {input.raw_invariant_vcf} \
		--intervals {params.chr} \
		-filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \
		-filter \"MQ < 40.0\" --filter-name \"MQ40\" \
		-filter \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" \
		-filter \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8\" \
		-G-filter \"DP < 6 || DP > 200\" \
		-G-filter-name \"DP_6-200\" \
		-O {output.invariant_flagged}")

rule filtered_to_nocall_invariant:
	input:
		invariant_flagged = "data/interim/filtering/{chr}_allsamps.filtered.qual.dp6_200.invariantOnly.vcf",
		invariant_idx = "data/interim/filtering/{chr}_allsamps.filtered.qual.dp6_200.invariantOnly.vcf.idx"
	output:
		invariant_nocall = temp("data/interim/filtering/{chr}_allsamps.filtered.qual.dp6_200.invariant.nocall.vcf"),
		invariant_nocall_idx = temp("data/interim/filtering/{chr}_allsamps.filtered.qual.dp6_200.invariant.nocall.vcf.idx")
	run:
		shell("gatk VariantFiltration \
		-V {input.invariant_flagged} \
		--set-filtered-genotype-to-no-call true \
		-O {output.invariant_nocall}")

rule filter_invariant:
	input:
		invariant_nocall = "data/interim/filtering/{chr}_allsamps.filtered.qual.dp6_200.invariant.nocall.vcf",
		invariant_nocall_idx = "data/interim/filtering/{chr}_allsamps.filtered.qual.dp6_200.invariant.nocall.vcf.idx",
		ref = config['ref']
	output:
		filtered_invariant = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp6_200.maxnocall10.invariant.sites.vcf",
		filtered_invariant_idx = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp6_200.maxnocall10.invariant.sites.vcf.idx"
	run:
		shell("gatk SelectVariants \
		-R {input.ref} \
		-V {input.invariant_nocall} \
		--max-nocall-fraction 0.1 \
		--exclude-filtered true \
		-O {output.filtered_invariant}")

# rule filter_snps:
# 	input:
# 		ref = config['ref'],
# 		raw_snp_vcf = "data/raw/vcf_bpres/all_samps.raw.snps.vcf",
# 	output:
# 		filtered_snps = pipe("data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp6_200.maxnocall10.biallelic.snps.vcf"),
# 		filtered_snps_indx = pipe("data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp6_200.maxnocall10.biallelic.snps.vcf.idx")
# 	params:
# 		chr = "{chr}",
# 		scratch_dir = config['scratch'],
# 		filter_tmp = config['scratch'] + "{chr}_allsamps.filtered.qual.dp6_200.snpsOnly.vcf",
# 		nocall_tmp = config['scratch'] + "{chr}_allsamps.filtered.qual.dp6_200.snpsOnly.nocall.vcf",
# 		nocall_frac_tmp = config['scratch'] + "data/interim/{chr}_allsamps.filtered.qual.dp6_200.maxnocall10.biallelic.snps.vcf"
# 	run:
# 		shell("mkdir -p {params.scratch_dir}")
# 		shell("gatk VariantFiltration \
# 		-V {input.raw_snp_vcf} \
# 		--intervals {params.chr} \
# 		-filter \"QD < 2.0\" --filter-name \"QD2\" \
# 		-filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \
# 		-filter \"SOR > 3.0\" --filter-name \"SOR3\" \
# 		-filter \"FS > 60.0\" --filter-name \"FS60\" \
# 		-filter \"MQ < 40.0\" --filter-name \"MQ40\" \
# 		-filter \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" \
# 		-filter \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8\" \
# 		-G-filter \"DP < 6 || DP > 200\" \
# 		-G-filter-name \"DP_6-200\" \
# 		-O {params.filter_tmp}")
# 		shell("gatk VariantFiltration \
# 		-V {params.filter_tmp} \
# 		--set-filtered-genotype-to-no-call true \
# 		-O {params.nocall_tmp}")
# 		shell("rm -rf {params.filter_tmp}*")
# 		shell("gatk SelectVariants \
# 		-R {input.ref} \
# 		-V {params.nocall_tmp} \
# 		-select-type SNP \
# 		--restrict-alleles-to BIALLELIC \
# 		--max-nocall-fraction 0.1 \
# 		--exclude-filtered true \
# 		-O {params.nocall_frac_tmp}")
# 		shell("rm -rf {params.nocall_tmp}*")
# 		shell("mv {params.nocall_frac_tmp} {output.filtered_snps}")
# 		shell("mv {params.nocall_frac_tmp}.idx {output.filtered_snps_indx}")
# 		shell("rm -rf {params.nocall_frac_tmp}*")


# rule filter_invariant:
# 	input:
# 		ref = config['ref'],
# 		invariant_vcf = "data/raw/vcf_bpres/all_samps.raw.invariant.vcf"
# 	output:
# 		filtered_invariant = temp("data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp6_200.maxnocall10.invariant.vcf"),
# 		filtered_invariant_indx = temp("data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp6_200.maxnocall10.invariant.vcf.idx")
# 	params:
# 		chr = "{chr}",
# 		scratch_dir = "data/interim/filtered_vcf_bpres",
# 		filter_tmp = "data/interim/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp6_200.invariantOnly.vcf",
# 		nocall_tmp = "data/interim/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp6_200.invariantOnly.nocall.vcf",
# 		nocall_frac_tmp = "data/iterim/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp6_200.maxnocall10.invariant.vcf"
# 	run:
# 		shell("mkdir -p {params.scratch_dir}")
# 		shell("gatk VariantFiltration \
# 		-V {input.invariant_vcf} \
# 		--intervals {params.chr} \
# 		-filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \
# 		-filter \"MQ < 40.0\" --filter-name \"MQ40\" \
# 		-filter \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" \
# 		-filter \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8\" \
# 		-G-filter \"DP < 6 || DP > 200\" \
# 		-G-filter-name \"DP_6-200\" \
# 		-O {params.filter_tmp}")
# 		shell("gatk VariantFiltration \
# 		-V {params.filter_tmp} \
# 		--set-filtered-genotype-to-no-call true \
# 		-O {params.nocall_tmp}")
# 		shell("rm -rf {params.filter_tmp}*")
# 		shell("gatk SelectVariants \
# 		-R {input.ref} \
# 		-V {params.nocall_tmp} \
# 		--max-nocall-fraction 0.1 \
# 		--exclude-filtered true \
# 		-O {params.nocall_frac_tmp}")
# 		shell("rm -rf {params.nocall_tmp}*")
# 		shell("mv {params.nocall_frac_tmp} {output.filtered_invariant}")
# 		shell("mv {params.nocall_frac_tmp}.idx {output.filtered_invariant_indx}")
# 		shell("rm -rf {params.nocall_frac_tmp}*")

rule merge_filtered_snps:
	input:
		snps = expand("data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp6_200.maxnocall10.biallelic.snps.vcf", chr = CHR)
	output:
		merged_snps = "data/processed/filtered_vcf_bpres/allsamps.filtered.qual.dp6_200.maxnocall10.biallelic.snps.vcf.gz"
	params:
		input = lambda wildcards, input: " -I ".join(input)
	run:
		shell("gatk MergeVcfs \
		-I={params.input} \
		-O={output.merged_snps}")

rule merge_filtered_vcf:
	input:
		snps = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp6_200.maxnocall10.biallelic.snps.vcf",
		invariant = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp6_200.maxnocall10.invariant.vcf"
	output:
		merged = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp6_200.maxnocall10.allsites.vcf.gz"
	params:
		scratch_dir = config['scratch'],
		tmp = config['scratch'] + "{chr}_allsamps.filtered.qual.dp6_200.maxnocall10.allsites.vcf.gz"
	run:
		shell("mkdir -p {params.scratch_dir}")
		shell("gatk MergeVcfs \
		-I={input.snps} \
		-I={input.invariant} \
		-O={params.tmp}")
		shell("mv {params.tmp} {output.merged}")
		shell("mv {params.tmp}.idx {output.merged}.idx")
		shell("rm -rf {params.tmp}*")

# filter for depth
# convert filtered genotypes to no call
# replaces flagged genotypes with ./.

# rule filter_nocall:
# 	input:
# 		qual_vcf = "data/processed/filtered_snps_bpres/{chr}_all_samps.filtered.snps.invariant.vcf.gz",
# 		ref = config['ref']
# 	output:
# 		filtered_vcf = "data/processed/filtered_snps_bpres/{chr}_all_samps.filtered.snps.invariant.dp6_200.maxnocall10.vcf.gz"
# 	params:
# 		tmp_vcf = temp(config['scratch'] + "filtered_snps_bpres/{chr}_all_samps.filtered.snps.invariant.dp6_200.nocall.vcf.gz")
# 	run:
# 		shell("gatk VariantFiltration \
# 		-V {input.qual_vcf} \
# 		-G-filter \"DP < 6 || DP > 200\" \
# 		-G-filter-name \"DP_6-200\" \
# 		--set-filtered-genotype-to-no-call true \
# 		-O {params.tmp_vcf}")
# 		shell("gatk SelectVariants \
# 		-R {input.ref} \
# 		-V {params.tmp_vcf} \
# 		-select-type SNP \
# 		--restrict-alleles-to BIALLELIC \
# 		--max-nocall-fraction 0.1 \
# 		--exclude-filtered true \
# 		-O {output.filtered_vcf}")
# 		shell("rm -rf {params.tmp_vcf}")
		# shell("gatk VariantFiltration \
		# -V {input} \
		# --set-filtered-genotype-to-no-call true \
		# -O {output}")

# rule snps_only:
# 	input:
# 		ref = config['ref'],
# 		vcf = "data/processed/filtered_snps_bpres/{chr}_all_samps.filtered.snps.invariant.dp6_200.maxnocall10.vcf.gz"
# 	output:
# 		"data/processed/filtered_snps_bpres/{chr}_all_samps.filtered.snpsOnly.dp6_200.maxnocall10.vcf.gz"
# 	run:
# 		shell("gatk SelectVariants \
# 		-R {input.ref} \
# 		-V {input.vcf} \
# 		-select-type SNP \
# 		--restrict-alleles-to BIALLELIC \
# 		--max-nocall-fraction 0.1 \
# 		--exclude-filtered true \
# 		-O {output}")

# # add maxnocall fraction in sequential filtering step
# rule filter_nocall:
# 	input:
# 		ref = "data/external/ref/Boleracea_chromosomes.fasta",
# 		vcf = "data/processed/filtered_snps_bpres/{chr}.filtered.snps.vcf"
# 	output:
# 		"data/processed/filtered_snps_bpres/{chr}.filtered.nocall.vcf"
# 	run:
# 		shell("gatk SelectVariants \
# 		-V {input.vcf} \
# 		--max-nocall-fraction 0 \
# 		--exclude-filtered true \
# 		--restrict-alleles-to BIALLELIC \
# 		-O {output}")
#
# rule filter_depth:
# 	input:
# 		ref = "data/external/ref/Boleracea_chromosomes.fasta",
# 		vcf = "data/processed/filtered_snps_bpres/{chr}.filtered.nocall.vcf"
# 	output:
# 		dp = "data/processed/filtered_snps_bpres/{chr}.filtered.dp6_200.snps.vcf",
# 		dp2 = "data/processed/filtered_snps_bpres/{chr}.filtered.dp6_200.nocall.snps.vcf"
# 	run:
# 		shell("gatk VariantFiltration \
# 		-V {input.vcf} \
# 		-G-filter \"DP < 6 || DP > 200\" \
# 		-G-filter-name \"DP_6-200\" \
# 		--set-filtered-genotype-to-no-call true \
# 		-O {output.dp}")
# 		shell("gatk SelectVariants \
# 		-V {output.dp} \
# 		--max-nocall-fraction 0.1 \
# 		--exclude-filtered true \
# 		--restrict-alleles-to BIALLELIC \
# 		-O {output.dp2}")
#
# rule depth:
# 	input:
# 		vcf = "data/processed/filtered_snps_bpres/{chr}.filtered.dp6_200.nocall.snps.vcf",
# 		ref = "data/external/ref/Boleracea_chromosomes.fasta"
# 	output:
# 		dp = "reports/filtering_bpres/dp_{chr}.filtered.dp6_200"
# 	run:
# 		shell("gatk VariantsToTable \
# 		-R {input.ref} \
# 		-V {input.vcf} \
# 		-F CHROM -F POS -GF DP \
# 		-O {output.dp}")

# rule bgzip_vcf:
#     input:
#         "data/processed/filtered_snps_bpres/{chr}.filtered.dp6_200.nocall.snps.vcf"
#     output:
#         "data/processed/filtered_snps_bpres/{chr}.filtered.dp6_200.nocall.snps.vcf.gz"
#     run:
#         shell("bgzip {input}")
#         shell("tabix -p vcf {output}")
#
# rule combine_vcfs:
# 	input:
# 		expand("data/processed/filtered_snps_bpres/{chr}.filtered.dp6_200.nocall.snps.vcf.gz", chr = chr)
# 	output:
# 		"data/processed/filtered_snps_bpres/oleracea_filtered.vcf.gz"
# 	run:
# 		shell("bcftools concat {input} -Oz -o {output}")

# evaluate filtering parameters
# https://software.broadinstitute.org/gatk/documentation/article.php?id=11069
# not sure this is possible for this data set... need dbSNP info


# depth of coverage
# rule coverage_depth:
#     input:
#         all = expand("data/raw/sorted_reads/{sample}.sorted.bam", sample = SAMPLES),
#         bamlist = "models/angsd/ALL.bamlist",
#         ref = "data/external/ref/Boleracea_chromosomes.fasta"
#     output:
#         "reports/coverage/{chr}.coverage.txt"
#     params:
#         chr = "{chr}"
#     shell:
#         "samtools depth \
#         -f models/angsd/ALL.bamlist \
#         -q 20 \
#         -Q 30 \
#         -r {params.chr} \
#         --reference {input.ref} > {output}"
