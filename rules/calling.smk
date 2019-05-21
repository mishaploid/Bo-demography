# Step 2: variant discovery!

# run GATK haplotype caller

rule hap_caller:
	input:
		ref = "data/external/ref/Boleracea_chromosomes.fasta",
		bam = "data/interim/add_rg/{sample}.rg.dedup.bam",
		bai = "data/interim/add_rg/{sample}.rg.dedup.bai"
	output:
		"data/interim/{sample}.raw.snps.indels.g.vcf"
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
		--emit-ref-confidence GVCF")

# replace sample names

# rule rename_samples:
# 	input:
# 		"data/interim/{sample}.raw.snps.indels.g.vcf"
# 	output:
# 		"data/interim/{sample}.renamed.raw.snps.indels.g.vcf"
# 	params:
# 		sample = "{sample}"
# 	run:
# 		shell("gatk RenameSampleInVcf \
# 		-I={input} \
# 		-O={output} \
# 		--NEW_SAMPLE_NAME={params.sample} \
# 		--CREATE_INDEX=true")

# combine GVCFs
# https://software.broadinstitute.org/gatk/documentation/article?id=11813
# snakemake considerations - https://bitbucket.org/snakemake/snakemake/issues/895/combine-multiple-files-for-input-but

rule combine_gvcfs:
	input:
		expand("data/interim/{sample}.raw.snps.indels.g.vcf", sample = SAMPLES)
	output:
		directory("data/interim/combined_database/{chr}")
	params:
		files = lambda wildcards, input: " -V ".join(input),
		region = "{chr}"
	run:
		shell("gatk GenomicsDBImport \
		-V {params.files} \
		--genomicsdb-workspace-path {output} \
		--batch-size 50 \
		--intervals {params.region}")

# joint genotyping - raw SNP and indel VCF

rule joint_geno:
	input:
		dir = "data/interim/combined_database/{chr}/vcfheader.vcf",
		ref = "data/external/ref/Boleracea_chromosomes.fasta"
	output:
		"data/raw/{chr}.raw.snps.indels.vcf"
	params:
		db = "gendb://data/interim/combined_database/{chr}"
	run:
		shell("gatk GenotypeGVCFs \
		-R {input.ref} \
		-V {params.db} \
		-new-qual \
		-G StandardAnnotation \
		-G AS_StandardAnnotation \
		-O {output}")
