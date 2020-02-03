################################################################################
## Rules for mapping short reads to a reference
## obtain reference/fasterq-dump -> trim reads (trimmomatic) -> align (bwa-mem)
## -> sort bam files (SortSam) -> mark PCR duplicates (MarkDuplicates)
################################################################################

################################################################################
## Get reference
# see https://gatkforums.broadinstitute.org/gatk/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk
# download reference (wget)
# generate BWA index (bwa index)
#    collection of files used by BWA for alignment
# generate fasta index (samtools faidx) -
#    .fai file with one record per line for each contig in FASTA reference
#    contig name | size | location | basesPerLine | bytesPerLine
# create sequence dictionary (gatk CreateSequenceDictionary)
#    .dict file formatted like SAM header that describes contents of FASTA reference
################################################################################

rule get_ref:
    output:
        ref = "data/external/ref/Boleracea_chromosomes.fasta",
        dict = "data/external/ref/Boleracea_chromosomes.dict"
    run:
        shell("wget --directory-prefix=data/external/ref \
        http://www.genoscope.cns.fr/externe/plants/data/Boleracea_chromosomes.fasta")
        shell("bwa index -a bwtsw {output.ref}")
        shell("samtools faidx {output.ref}")
        shell("gatk CreateSequenceDictionary \
        -R={output.ref} \
        -O={output.dict}")

################################################################################
# OPTION 1: FASTQ to BAM all in one step
# Download FASTQ files (fasterq-dump)
#   uses a temp directory
# run Trimmomatic (consider removing this step)
#   removes leading/trailing low quality bases (LEADING/TRAILING)
#   scans with a 4-base sliding window and cuts when quality per base < 15 (SLIDINGWINDOW)
#   removes reads less than 36 bases long (MINLEN)
# align with bwa-mem
################################################################################

rule fastq2bam:
    input:
        samp_ids = "data/external/Sra_oleracea.csv",
        ref = "data/external/ref/Boleracea_chromosomes.fasta"
    output:
        temp(touch("data/interim/mapped_reads/{sample}.bam"))
    params:
        tmp = "/scratch/sdturner/map_reads/{sample}",
        sample = "{sample}",
        SRR = lambda wildcards: sample_dict[wildcards.sample],
        stem = "/scratch/sdturner/map_reads/{sample}/{sample}"
    wildcard_constraints:
        sample = "B_cretica_[A-D]"
    threads: 24
    run:
        print({params.sample}, {params.SRR})
        shell("mkdir -p {params.tmp}")
        shell("fasterq-dump \
        {params.SRR} \
        -O {params.tmp} \
        -o {params.sample} \
        -t {params.tmp} \
        -e {threads} \
        -p")
        shell("java -jar /share/apps/Trimmomatic-0.36/trimmomatic.jar PE \
        {params.stem}_1.fastq {params.stem}_2.fastq \
        {params.stem}.forward.1.fastq \
        {params.stem}.unpaired.1.fastq \
        {params.stem}.reverse.2.fastq \
        {params.stem}.unpaired.2.fastq \
        LEADING:3 \
		  TRAILING:3 \
		  SLIDINGWINDOW:4:15 \
		  MINLEN:36")
        shell("(bwa mem -t {threads} \
        {input.ref} \
        {params.stem}.forward.1.fastq \
        {params.stem}.reverse.2.fastq | \
        samtools view -Sb > {output})")
        shell("rm -rf {params.tmp}")

################################################################################
# OPTION 2: FASTQ files stored locally
# run Trimmomatic (consider removing this step)
#   removes leading/trailing low quality bases (LEADING/TRAILING)
#   scans with a 4-base sliding window and cuts when quality per base < 15 (SLIDINGWINDOW)
#   removes reads less than 36 bases long (MINLEN)
# align with bwa-mem
################################################################################

# TODO: add readgroup information at this step using tip from Michelle:
# bwa mem -R $(echo "@RG\tID:${LANE}_${SAMPLE}\tSM:$SAMPLE\tLB:${SAMPLE}.1\tPL:ILLUMINA") -a $B73FA $R1 $R2

rule bwa_mem:
    input:
        fastq1 = "data/raw/fastq/{sample}_R1_001.fastq.gz",
        fastq2 = "data/raw/fastq/{sample}_R2_001.fastq.gz",
        ref = "data/external/ref/Boleracea_chromosomes.fasta"
    output:
        temp(touch("data/interim/mapped_reads/{sample}.bam"))
    threads: 6
    params:
        tmp = "/scratch/sdturner/map_reads/{sample}",
        stem = "/scratch/sdturner/map_reads/{sample}/{sample}"
    wildcard_constraints:
        pattern = "!B_cretica_[A-D]"
    run:
        shell("mkdir -p {params.tmp}")
        shell("java -jar /share/apps/Trimmomatic-0.36/trimmomatic.jar PE \
        {input.fastq1} {input.fastq2} \
        {params.stem}.forward.1.fastq \
        {params.stem}.unpaired.1.fastq \
        {params.stem}.reverse.2.fastq \
        {params.stem}.unpaired.2.fastq \
		LEADING:3 \
		TRAILING:3 \
		SLIDINGWINDOW:4:15 \
		MINLEN:36")
        shell("(bwa mem \
        -t {threads} \
        {input.ref} \
        {params.stem}.forward.1.fastq \
        {params.stem}.reverse.2.fastq | \
        samtools view -Sb > {output})")
        shell("rm -rf {params.tmp}")

################################################################################
# sort BAM files with Picard
# uses a temp directory to speed things up
# creates an index file for each bam (--CREATE_INDEX=true)
################################################################################

rule sort_bam:
    input:
        "data/interim/mapped_reads/{sample}.bam"
    output:
    	"data/raw/sorted_reads/{sample}.sorted.bam"
    params:
    	tmp = "/scratch/sdturner/sort_bam/{sample}"
    run:
        shell("mkdir -p {params.tmp}")
    	shell("gatk SortSam \
    	-I={input} \
    	-O={output} \
    	--SORT_ORDER=coordinate \
    	--TMP_DIR={params.tmp} \
    	--CREATE_INDEX=true")
        shell("rm -rf {params.tmp}")

################################################################################
# mark duplicates with Picard
# note: no need to remove duplicates here - Haplotype Caller will ignore them
# returns qc metrics on % optical duplicates & histogram of return on coverage
#   e.g. if you were to sequence BIN times more than current level,
#   you would see VALUE times more coverage
# https://github.com/broadinstitute/picard/issues/917
# http://broadinstitute.github.io/picard/faq.html
# duplicate reads are flagged with 0x0400
# could also remove with 'samtools view -F 0x400 sorted.bam'
# metrics file - records duplication metrics (required)
# create_index - true; indexes resulting BAM file
# -MAX_FILE_HANDLES - reduced this to meet system requirements for max number of open files (check with `ulimit -n`)
################################################################################

rule mark_dups:
    input:
        "data/raw/sorted_reads/{sample}.sorted.bam"
    output:
        bam = temp(touch("data/interim/mark_dups/{sample}.dedup.bam")),
        index = temp(touch("data/interim/mark_dups/{sample}.dedup.bai")),
        metrics = "qc/mark_dup/{sample}_metrics.txt"
    params:
        tmp = "/scratch/sdturner/mark_dups/{sample}"
    run:
        shell("mkdir -p {params.tmp}")
        shell("gatk MarkDuplicates \
        -I={input} \
        -O={output.bam} \
        --METRICS_FILE={output.metrics} \
        --CREATE_INDEX=true \
        -MAX_FILE_HANDLES=1000 \
        --ASSUME_SORT_ORDER=coordinate \
        --TMP_DIR={params.tmp}")
        shell("rm -rf {params.tmp}")

################################################################################
# add read groups with Picard
# RGID, RGLB, RGPL, and RGPU are dummy variables
# RGSM adds the sample id to each bam file
################################################################################
rule add_rg:
    input:
    	"data/interim/mark_dups/{sample}.dedup.bam"
    output:
    	bam = temp(touch("data/interim/add_rg/{sample}.rg.dedup.bam")),
    	index = temp(touch("data/interim/add_rg/{sample}.rg.dedup.bai"))
    params:
        sample = "{sample}",
        tmp = "/scratch/sdturner/add_rg/{sample}"
    run:
        shell("mkdir -p {params.tmp}")
    	shell("gatk AddOrReplaceReadGroups \
    	-I={input} \
    	-O={output.bam} \
    	--CREATE_INDEX=true \
    	-RGID=4 \
    	-RGLB=lib1 \
    	-RGPL=illumina \
    	-RGPU=unit1 \
    	-RGSM={params.sample} \
    	--TMP_DIR {params.tmp}")
        shell("rm -rf {params.tmp}")

################################################################################
# quality metrics with qualimap
# don't freak out about insert size distribution (it's a bwa-mem thing):
# https://github.com/lh3/bwa/issues/113
################################################################################

rule bamqc:
	input:
		bam = "data/raw/sorted_reads/{sample}.sorted.bam",
		gff = "data/external/ref/Boleracea_annotation.gff"
	output:
		"reports/bamqc/{sample}_stats/qualimapReport.html"
	params:
		dir = "reports/bamqc/{sample}_stats"
	threads: 8
	run:
		shell("qualimap bamqc \
		-bam {input.bam} \
		--paint-chromosome-limits \
		-gff {input.gff} \
		-nt {threads} \
		-outdir {params.dir} \
		-outformat HTML \
		--skip-duplicated")

################################################################################
# Combine qualimap results
# this includes some python magic to create the params.infile txt file, which
# lists all of the individual qualimap report information
################################################################################

rule multibamqc:
	input:
		all = expand("reports/bamqc/{sample}_stats/qualimapReport.html", sample = SAMPLES)
	output:
		"reports/multisampleBamQcReport.html"
	params:
		outdir = "reports",
		infile = "models/bamqc_list.txt"
	run:
		shell("find reports/bamqc -mindepth 1 -maxdepth 1 -type d | grep SamC > models/ALL.bamqclist.txt")
		import pandas as pd
		data = pd.read_csv("models/ALL.bamqclist.txt", sep = " ", header = None, names = ['filename'])
		data['sample'] = data['filename'].str.split('.').str[0].str.split('/').str[2].str.split('_stats').str[0]
		data = data.sort_values('sample', axis = 0, ascending = True)
		data = data[['sample','filename']]
		data.to_csv(r'models/bamqc_list.txt', header = None, index = None, sep = ' ', mode = 'a')
		shell("qualimap multi-bamqc \
		--data {params.infile} \
		--paint-chromosome-limits \
		-outdir {params.outdir} \
		-outformat html \
		--java-mem-size=40G")
