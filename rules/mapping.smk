### Data cleanup

# trim adapters
# rule trimmomatic_pe:
#     input:
#         temp("data/external/fastq_raw/{sample}_pass_1.fastq.gz"),
#         temp("data/external/fastq_raw/{sample}_pass_2.fastq.gz")
#     output:
#         temp("data/interim/trimmed_reads/{sample}.1.fastq.gz"),
#         temp("data/interim/trimmed_reads/{sample}.1.unpaired.fastq.gz"),
#         temp("data/interim/trimmed_reads/{sample}.2.fastq.gz"),
#         temp("data/interim/trimmed_reads/{sample}.2.unpaired.fastq.gz")
#     log:
#         "logs/trimmomatic/{sample}.log"
#     run:
#         shell("java -jar {trimmomatic} PE {input} {output} \
# 		LEADING:3 \
# 		TRAILING:3 \
# 		SLIDINGWINDOW:4:15 \
# 		MINLEN:36 2> {log}")

# LEADING:3 - remove leading low quality or N bases (quality < 3)
# TRAILING:3 - remove trailing low quality or N bases (quality < 3)
# SLIDINGWINDOW: scan read with 4-base sliding window, cutting when average quality
# drops below 15
# MINLEN:36 - drop reads below 36 bases long
# note: check sizes of paired vs. unpaired output files (expect paired to be larger)

# map to reference (TO1000 v2.1 from Parkin et al. 2014)
# rule bwa_map:
#     input:
#        "data/external/ref/Brassica_oleracea.v2.1.dna.toplevel.fa",
#        temp("data/interim/trimmed_reads/{sample}.1.fastq.gz"),
#        temp("data/interim/trimmed_reads/{sample}.2.fastq.gz")
#     output:
#        temp("data/interim/mapped_reads/{sample}.bam")
# 	log:
#        "logs/bwa/{sample}.log"
# 	threads: 8
# 	run:
# 		shell("(bwa mem -t {threads} {input} | \
# 		samtools view -Sb > {output}) 2> {log}")

# sort BAM files with Picard
rule sort_bam:
	input:
		"data/interim/mapped_reads/{sample}.bam"
	output:
		bam = "data/raw/sorted_reads/{sample}.sorted.bam",
		tmp = temp(directory("data/sort_bam/{sample}"))
	run:
		shell("gatk SortSam \
		-I={input} \
		-O={output.bam} \
		--SORT_ORDER=coordinate \
		--TMP_DIR={output.tmp}")

# mark duplicates with Picard
# no need to remove duplicates here - Haplotype Caller will ignore them
rule mark_dups:
    input:
    	"data/raw/sorted_reads/{sample}.sorted.bam"
    output:
        bam = temp("data/raw/{sample}.dedup.bam"),
        index = temp("data/raw/{sample}.dedup.bai"),
        metrics = "qc/dedup_reads/{sample}_metrics.txt",
        tmp = temp(directory("data/mark_dups/{sample}"))
    run:
        shell("gatk MarkDuplicates \
        -I={input} \
        -O={output.bam} \
        --METRICS_FILE={output.metrics} \
        --CREATE_INDEX=true \
        -MAX_FILE_HANDLES=1000 \
        --ASSUME_SORT_ORDER=coordinate \
        --TMP_DIR={output.tmp}")

# metrics file - records duplication metrics (required)
# create_index - true; indexes resulting BAM file
# -MAX_FILE_HANDLES - reduced this to meet system requirements for max number of open files (check with `ulimit -n`)

# add read groups with Picard
# info is just dummy variables
rule add_rg:
	input:
		"data/raw/{sample}.dedup.bam"
	output:
		bam = temp("data/raw/{sample}.rg.dedup.bam"),
		index = temp("data/raw/{sample}.rg.dedup.bai"),
		tmp = temp(directory("data/add_rg/{sample}"))
	run:
		shell("gatk AddOrReplaceReadGroups \
		-I={input} \
		-O={output.bam} \
		--CREATE_INDEX=true \
		-RGID=4 \
		-RGLB=lib1 \
		-RGPL=illumina \
		-RGPU=unit1 \
		-RGSM=20 \
		--TMP_DIR {output.tmp}")
