# hard filtering
# https://software.broadinstitute.org/gatk/documentation/article.php?id=2806
# Note: Variant Quality Score Recalibration is not an option (need truth/training sets)

# extract SNPs

rule get_snps:
	input:
		ref = "data/external/ref/Boleracea_chromosomes.fasta",
		vcf = "data/raw/{chr}.raw.snps.indels.vcf"
	output:
		"data/raw/{chr}.raw.snps.vcf"
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
		vcf = "data/raw/{chr}.raw.snps.vcf"
	output:
		"data/processed/{chr}.filtered.snps.vcf"
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
		-O {output}")

# evaluate filtering parameters
# https://software.broadinstitute.org/gatk/documentation/article.php?id=11069
# not sure this is possible for this data set... need dbSNP info

# quality metrics with qualimap
rule bamqc:
	input:
		all = expand("data/raw/sorted_reads/{sample}.sorted.bam", sample = SAMPLES)
	output:
		"reports/bamqc/report.pdf"
	params:
		outdir = "reports/bamqc",
		infile = "models/bamqc_list.txt"
	run:
		shell("ls -rt data/raw/sorted_reads/*.bam > models/ALL.bamlist.txt")
		import pandas as pd
		data = pd.read_csv("models/ALL.bamlist.txt", sep = " ", header = None, names = ['filename'])
		data['sample'] = data['filename'].str.split('.').str[0].str.split('/').str[3]
		data = data.drop([119], axis = 0)
		data = data.sort_values('sample', axis = 0, ascending = True)
		data = data[['sample','filename']]
		data.to_csv(r'models/bamqc_list.txt', header = None, index = None, sep = ' ', mode = 'a')
		shell("qualimap multi-bamqc \
		--data {params.infile} \
		--paint-chromosome-limits \
		-outdir {params.outdir} \
		-outformat PDF \
		--run-bamqc \
		--java-mem-size=60G")

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
