# hard filtering
# https://software.broadinstitute.org/gatk/documentation/article.php?id=2806
# Note: Variant Quality Score Recalibration is not an option (need truth/training sets)

# extract SNPs

rule get_snps:
	input:
		ref = "data/external/ref/Boleracea_chromosomes.fasta",
		vcf = "data/raw/vcf/{chr}.raw.snps.indels.vcf"
	output:
		"data/raw/vcf/{chr}.raw.snps.vcf"
	run:
		shell("gatk SelectVariants \
		-R {input.ref} \
		-V {input.vcf} \
		-select-type SNP \
		-O {output}")

# apply filters
# using base GATK recommendations
# https://software.broadinstitute.org/gatk/documentation/article?id=23216#2

rule filter_snps:
	input:
		ref = "data/external/ref/Boleracea_chromosomes.fasta",
		vcf = ancient("data/raw/vcf/{chr}.raw.snps.vcf")
	output:
		"data/processed/filtered_snps/{chr}.filtered.snps.vcf"
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

rule bgzip_vcf:
    input:
        "data/processed/filtered_snps/{chr}.filtered.snps.vcf"
    output:
        "data/processed/filtered_snps/{chr}.filtered.snps.vcf.gz"
    run:
        shell("bgzip {input}")
        shell("tabix -p vcf {output}")

rule combine_vcfs:
	input:
		expand("data/processed/filtered_snps/{chr}.filtered.snps.vcf.gz", chr = chr)
	output:
		"data/processed/filtered_snps/oleracea_combined.vcf.gz"
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
