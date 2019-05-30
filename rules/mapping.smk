# Rules for mapping
# obtain reference/fasterq-dump -> trim reads (trimmomatic) -> align (bwa-mem)
# sort bam files (SortSam)
# mark PCR duplicates (MarkDuplicates)

# download reference (wget)
# generate BWA index (bwa index)
# generate fasta index (samtools faidx)
# generate sequence dictionary (gatk CreateSequenceDictionary)
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

# use fastq dump to download FASTQ files
# get SRA info for B. oleracea:
# cat SraRunInfo.csv | grep oleracea > Sra_oleracea.csv
# sed -i "1s/^/$(head -n1 SraRunInfo.csv)\n/" Sra_oleracea.csv

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

# sort BAM files with Picard
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

# mark duplicates with Picard
# no need to remove duplicates here - Haplotype Caller will ignore them
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

# metrics file - records duplication metrics (required)
# create_index - true; indexes resulting BAM file
# -MAX_FILE_HANDLES - reduced this to meet system requirements for max number of open files (check with `ulimit -n`)

# add read groups with Picard
# info is just dummy variables
### TODO: test if changing -RGSM flag fixes sample naming
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
