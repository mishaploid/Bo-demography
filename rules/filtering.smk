# hard filtering
# https://software.broadinstitute.org/gatk/documentation/article.php?id=2806
# Note: Variant Quality Score Recalibration is not an option (need truth/training sets)

# extract SNPs

rule get_snps:
	input:
		ref = "data/external/ref/Boleracea_chromosomes.fasta",
		vcf = "data/raw/vcf_bpres/{chr}.raw.snps.indels.vcf"
	output:
		"data/raw/vcf_bpres/{chr}.raw.snps.vcf"
	run:
		shell("gatk SelectVariants \
		-R {input.ref} \
		-V {input.vcf} \
		-select-type SNP \
		-O {output}")

# apply filters
# using base GATK recommendations
# https://software.broadinstitute.org/gatk/documentation/article?id=23216#2
# https://software.broadinstitute.org/gatk/documentation/article.php?id=3225

rule filter_snps:
	input:
		ref = "data/external/ref/Boleracea_chromosomes.fasta",
		vcf = ancient("data/raw/vcf_bpres/{chr}.raw.snps.vcf")
	output:
		"data/processed/filtered_snps_bpres/{chr}.filtered.snps.vcf"
	run:
		shell("gatk VariantFiltration \
		-V {input.vcf} \
		-filter \"QD < 2.0\" --filter-name \"QD2\" \
		-filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \
		-filter \"SOR > 3.0\" --filter-name \"SOR3\" \
		-filter \"FS > 60.0\" --filter-name \"FS60\" \
		-filter \"MQ < 40.0\" --filter-name \"MQ40\" \
		-filter \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" \
		-filter \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8\" \
		--set-filtered-genotype-to-no-call true \
		-O {output}")

# filtering diagnostics
# https://evodify.com/gatk-in-non-model-organism/

rule diagnostics:
	input:
		vcf = "data/processed/filtered_snps_bpres/{chr}.filtered.snps.vcf",
		ref = "data/external/ref/Boleracea_chromosomes.fasta"
	output:
		# stats = "reports/filtering/gvcf_{chr}.table",
		filtered = "reports/filtering_bpres/gvcf_{chr}.filtered"
	run:
		shell("gatk VariantsToTable \
		-R {input.ref} \
		-V {input.vcf} \
		-F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F AD \
		-O {output.filtered}")
		# shell("gatk VariantsToTable \
		# -R {input.ref} \
		# -V {input.vcf} \
		# -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F AD \
		# -O {output.stats}")

# add maxnocall fraction in sequential filtering step
rule filter_nocall:
	input:
		ref = "data/external/ref/Boleracea_chromosomes.fasta",
		vcf = "data/processed/filtered_snps_bpres/{chr}.filtered.snps.vcf"
	output:
		"data/processed/filtered_snps_bpres/{chr}.filtered.nocall.vcf"
	run:
		shell("gatk SelectVariants \
		-V {input.vcf} \
		--max-nocall-fraction 0 \
		--exclude-filtered true \
		--restrict-alleles-to BIALLELIC \
		-O {output}")

rule filter_depth:
	input:
		ref = "data/external/ref/Boleracea_chromosomes.fasta",
		vcf = "data/processed/filtered_snps_bpres/{chr}.filtered.nocall.vcf"
	output:
		dp = "data/processed/filtered_snps_bpres/{chr}.filtered.dp6_200.snps.vcf",
		dp2 = "data/processed/filtered_snps_bpres/{chr}.filtered.dp6_200.nocall.snps.vcf"
	run:
		shell("gatk VariantFiltration \
		-V {input.vcf} \
		-G-filter \"DP < 6 || DP > 200\" \
		-G-filter-name \"DP_6-200\" \
		--set-filtered-genotype-to-no-call true \
		-O {output.dp}")
		shell("gatk SelectVariants \
		-V {output.dp} \
		--max-nocall-fraction 0.1 \
		--exclude-filtered true \
		--restrict-alleles-to BIALLELIC \
		-O {output.dp2}")

rule depth:
	input:
		vcf = "data/processed/filtered_snps_bpres/{chr}.filtered.dp6_200.nocall.snps.vcf",
		ref = "data/external/ref/Boleracea_chromosomes.fasta"
	output:
		dp = "reports/filtering_bpres/dp_{chr}.filtered.dp6_200"
	run:
		shell("gatk VariantsToTable \
		-R {input.ref} \
		-V {input.vcf} \
		-F CHROM -F POS -GF DP \
		-O {output.dp}")

rule bgzip_vcf:
    input:
        "data/processed/filtered_snps_bpres/{chr}.filtered.dp6_200.nocall.snps.vcf"
    output:
        "data/processed/filtered_snps_bpres/{chr}.filtered.dp6_200.nocall.snps.vcf.gz"
    run:
        shell("bgzip {input}")
        shell("tabix -p vcf {output}")

rule combine_vcfs:
	input:
		expand("data/processed/filtered_snps_bpres/{chr}.filtered.dp6_200.nocall.snps.vcf.gz", chr = chr)
	output:
		"data/processed/filtered_snps_bpres/oleracea_filtered.vcf.gz"
	run:
		shell("bcftools concat {input} -Oz -o {output}")

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
