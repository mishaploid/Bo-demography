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
