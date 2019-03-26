# Step 2: variant discovery!

# haplotype caller

rule hap_caller:
	input:
		ref = "data/external/ref/Brassica_oleracea.v2.1.dna.toplevel.fa",
		bam = "data/raw/{sample}.rg.dedup.bam",
		bai = "data/raw/{sample}.rg.dedup.bai"
	output:
		temp("data/sdturner/{sample}.raw.snps.indels.g.vcf")
	run:
		shell("gatk HaplotypeCaller \
		-I {input.bam} \
		-O {output} \
		-R {input.ref} \
		-G StandardAnnotation \
		-G AS_StandardAnnotation \
		--emit-ref-confidence GVCF")

# combine GVCFs
# https://software.broadinstitute.org/gatk/documentation/article?id=11813
# snakemake considerations - https://bitbucket.org/snakemake/snakemake/issues/895/combine-multiple-files-for-input-but
# expand("data/interim/{sample}.raw.snps.indels.g.vcf", sample = SAMPLES.remove("SRR7881031"))
rule combine_gvcfs:
	input:
		expand("data/sdturner/{sample}.raw.snps.indels.g.vcf", sample = SAMPLES)
	output:
		temp(directory("data/interim/combined_database"))
	params:
		files = lambda wildcards, input: " -V ".join(input)
	run:
		shell("gatk GenomicsDBImport \
		-V {params.files} \
		-G StandardAnnotation \
		-G AS_StandardAnnotation \
		--genomicsdb-workspace-path {output} \
		--intervals C1,C2,C3,C4,C5,C6,C7,C8,C9")

# joint genotyping - raw SNP and indel VCF
rule joint_geno:
	input:
		dir = directory("data/interim/combined_database"),
		ref = "data/external/ref/Brassica_oleracea.v2.1.dna.toplevel.fa"
	output:
		"data/raw/raw_variants.vcf"
	params:
		db = "gendb://data/interim/combined_database"
	threads: 10
	run:
		shell("gatk GenotypeGVCFs \
		-R {input.ref} \
		-V {params.db} \
		-nt {threads} \
		-newQual \
		-G StandardAnnotation \
		-G AS_StandardAnnotation \
		-O {output}")
