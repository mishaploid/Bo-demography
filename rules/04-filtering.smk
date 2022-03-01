################################################################################
## Rules to perform hard filtering for variant and invariant sites
## Outputs:
##		Merged VCF of jointly called SNPs
##		Merged VCF with variant + invariant sites (for nucleotide diversity)
## https://software.broadinstitute.org/gatk/documentation/article.php?id=2806
## Note: Variant Quality Score Recalibration is not an option
## 		(need truth/training sets)
## Follows this helpful guide on filtering for non-model organisms:
## 		https://evodify.com/gatk-in-non-model-organism/
################################################################################

################################################################################
## Select variant sites and restrict alleles to be biallelic
################################################################################

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

################################################################################
## Select invariant sites
################################################################################

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

################################################################################
## Plot filtering diagnostics (https://evodify.com/gatk-in-non-model-organism/)
## NOTE: some annotations cannot be calculated for invariant sites, because
## 		they require a mix of ref and alt reads(and no GQ reported)
##		https://gatkforums.broadinstitute.org/gatk/discussion/4711/genotype-qualities-qual-scores-for-all-sites
################################################################################

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

################################################################################
## Annotate variants with filtering criteria using base GATK recommendations
## Resources:
## https://software.broadinstitute.org/gatk/documentation/article?id=23216#2
## https://software.broadinstitute.org/gatk/documentation/article.php?id=3225
## https://informatics.fas.harvard.edu/whole-genome-resquencing-for-population-genomics-fastq-to-vcf.html#filtering
## https://gatk.broadinstitute.org/hc/en-us/articles/360035531012--How-to-Filter-on-genotype-using-VariantFiltration
################################################################################

## Filters:
## 	QD = quality by depth
## 		variant quality score divided by depth of alternate allele
## 	QUAL = SNP quality
## 		phred-scaled quality score; measure in confidence that there is a variation at a site
## 	SOR = StrandOddsRatio
## 		high values suggest strand bias
## 	FS = Fisher strand
## 		phred-scaled p-values for strand bias; higher values more likely to be false +
## 	MQ = RMS Mapping quality
## 		root mean square of mapping quality of reads across samples
## 	MQRankSum = Mapping quality rank sum test
## 		U-based z-approx from Mann-Whitney Rank Sum Test for MQ
## 	ReadPosRankSum
## 		U-based z-approx from Mann-Whitney Rank Sum Test for distance of alt allele from end of reads
## 		bases near end of read more likely to have errors
## 		if all reads with the alt allele are near the end may be errors

rule filter_flags_snps:
	input:
		raw_snp_vcf = "data/raw/vcf_bpres/all_samps.raw.snps.vcf"
	output:
		snps_flagged = temp("data/interim/filtering/{chr}_allsamps.filtered.qual.dp5_200.snpsOnly.vcf"),
		snps_idx = temp("data/interim/filtering/{chr}_allsamps.filtered.qual.dp5_200.snpsOnly.vcf.idx")
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
		-G-filter \"DP < 5 || DP > 200\" \
		-G-filter-name \"DP_5-200\" \
		-O {output.snps_flagged}")

rule filtered_to_nocall_snps:
	input:
		snps_flagged = "data/interim/filtering/{chr}_allsamps.filtered.qual.dp5_200.snpsOnly.vcf",
		snps_idx = "data/interim/filtering/{chr}_allsamps.filtered.qual.dp5_200.snpsOnly.vcf.idx"
	output:
		snps_nocall = temp("data/interim/filtering/{chr}_allsamps.filtered.qual.dp5_200.snpsOnly.nocall.vcf"),
		snps_nocall_idx = temp("data/interim/filtering/{chr}_allsamps.filtered.qual.dp5_200.snpsOnly.nocall.vcf.idx")
	run:
		shell("gatk VariantFiltration \
		-V {input.snps_flagged} \
		--set-filtered-genotype-to-no-call true \
		-O {output.snps_nocall}")

rule filter_snps:
	input:
		snps_nocall = "data/interim/filtering/{chr}_allsamps.filtered.qual.dp5_200.snpsOnly.nocall.vcf",
		snps_nocall_idx = "data/interim/filtering/{chr}_allsamps.filtered.qual.dp5_200.snpsOnly.nocall.vcf.idx",
		ref = config['ref']
	output:
		filtered_snps = temp("data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf"),
		filtered_snps_idx = temp("data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.idx")
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
		invariant_flagged = temp("data/interim/filtering/{chr}_allsamps.filtered.qual.dp5_200.invariantOnly.vcf"),
		invariant_idx = temp("data/interim/filtering/{chr}_allsamps.filtered.qual.dp5_200.invariantOnly.vcf.idx")
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
		-G-filter \"DP < 5 || DP > 200\" \
		-G-filter-name \"DP_5-200\" \
		-O {output.invariant_flagged}")

rule filtered_to_nocall_invariant:
	input:
		invariant_flagged = "data/interim/filtering/{chr}_allsamps.filtered.qual.dp5_200.invariantOnly.vcf",
		invariant_idx = "data/interim/filtering/{chr}_allsamps.filtered.qual.dp5_200.invariantOnly.vcf.idx"
	output:
		invariant_nocall = temp("data/interim/filtering/{chr}_allsamps.filtered.qual.dp5_200.invariant.nocall.vcf"),
		invariant_nocall_idx = temp("data/interim/filtering/{chr}_allsamps.filtered.qual.dp5_200.invariant.nocall.vcf.idx")
	run:
		shell("gatk VariantFiltration \
		-V {input.invariant_flagged} \
		--set-filtered-genotype-to-no-call true \
		-O {output.invariant_nocall}")

rule filter_invariant:
	input:
		invariant_nocall = "data/interim/filtering/{chr}_allsamps.filtered.qual.dp5_200.invariant.nocall.vcf",
		invariant_nocall_idx = "data/interim/filtering/{chr}_allsamps.filtered.qual.dp5_200.invariant.nocall.vcf.idx",
		ref = config['ref']
	output:
		filtered_invariant = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.invariant.sites.vcf",
		filtered_invariant_idx = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.invariant.sites.vcf.idx"
	run:
		shell("gatk SelectVariants \
		-R {input.ref} \
		-V {input.invariant_nocall} \
		--max-nocall-fraction 0.1 \
		--exclude-filtered true \
		-O {output.filtered_invariant}")

rule bgzip_vcf:
	input:
		vcf = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf",
		index = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.idx"
	output:
		"data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz",
		"data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz.tbi"
	run:
		shell("bgzip {input}")
		shell("tabix -p vcf {input}.gz")

rule merge_filtered_snps:
	input:
		snps = expand("data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz", chr = CHR)
	output:
		merged_snps = "data/processed/filtered_vcf_bpres/allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz"
	params:
		input = lambda wildcards, input: " -I ".join(input)
	run:
		shell("gatk MergeVcfs \
		-I={params.input} \
		-O={output.merged_snps}")

rule merge_filtered_allsites:
	input:
		snps = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz",
		invariant = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.invariant.sites.vcf"
	output:
		allsites_vcf = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.allsites.vcf.gz"
	run:
		shell("gatk MergeVcfs \
		-I={input.snps} \
		-I={input.invariant} \
		-O={output.allsites_vcf}")

rule ind_stats:
	input:
		merged_snps = "data/processed/filtered_vcf_bpres/allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz"
	output:
		missingness = "reports/filtered.qual.dp5_200.maxnocall10.imiss",
		heterozygosity = "reports/filtered.qual.dp5_200.maxnocall10.het"
	run:
		shell("vcftools \
		--gzvcf {input.merged_snps} \
		--missing-indv \
		--out reports/filtered.qual.dp5_200.maxnocall10")
		shell("vcftools \
		--gzvcf {input.merged_snps} \
		--het \
		--out reports/filtered.qual.dp5_200.maxnocall10")
