### msmc analysis

## generate mappability mask
# masks repetitive/poorly mapped regions (e.g. centromeres)

# SNPable program
# generates a mask for single-end reads of length 'k' and stringency 'r'
# default is k = 35
# Thanks to Li Wang for helpful code!
# Code also adapted from https://github.com/JackyHess/MSMC_analysis

# split reference into overlapping kmer sequences
# the first step splits the reference fasta into read lengths of 100 bp
# the pipe combines -l reads into a short read file
# aligns reads to the genome using bwa aln
# generates sam file for single end reads using bwa samse
# perl script from SNPabl to generate the raw (fasta) mask (gen_raw_mask.pl)
# use gen_mask function from SNPable to apply stringency (-r)
# use makeMappabilityMask.py to convert fasta to bed
## requires download of makeMappabilityMask.py from https://github.com/stschiff/msmc-tools
### edit input/output files

# how are results for paired end reads affected by assuming single reads for mask?
# should be more stringent?
# better to use sampe?

# note had to update gzip.open in makeMappabilityMask.py from 'w' to 'wt'

rule mappability_mask:
    input:
        "data/external/ref/Boleracea_chromosomes.fasta"
    output:
        "data/processed/mappability_masks/Boleracea_chr{chr}.mask.bed.gz"
    params:
        tmp = "/scratch/sdturner/{chr}",
        chr = "{chr}",
        newfasta = "data/external/ref/Boleracea_{chr}.fa",
        famask = "data/interim/masks/{chr}_mask_100_50.fa"
    threads: 16
    run:
        shell("mkdir -p {params.tmp}")
        shell("samtools faidx {input} {params.chr} > {params.newfasta}")
        shell("{snpable}/splitfa {params.newfasta} 100 | split -l 80000000 - {params.tmp}/x")
        shell("cat {params.tmp}/x* > {params.tmp}/simu_reads.fq")
        shell("bwa aln -R 1000000 -O 3 -E 3 -t {threads} {input} {params.tmp}/simu_reads.fq > {params.tmp}/simu_reads.sai")
        shell("bwa samse {input} {params.tmp}/simu_reads.sai {params.tmp}/simu_reads.fq > {params.tmp}/simu_reads.sam")
        shell("perl {snpable}/gen_raw_mask.pl {params.tmp}/simu_reads.sam > {params.tmp}/raw_mask_100.fa")
        shell("{snpable}/gen_mask -l 100 -r 0.5  {params.tmp}/raw_mask_100.fa > {params.famask}")
        shell("rm -rf {params.tmp}")
        shell("./src/makeMappabilityMask.py {params.famask}")

### individual masks for to account for coverage by individual
# estimates mean coverage for each individual
# uses msmc-tools helper script bamCaller.py to create a mask for each individual
# nested shell command reports depth of coverage
# note: no way to get consensus quality (FQ) in GATK vcf (required for bamCaller.py)
##  see thread here: https://sourceforge.net/p/samtools/mailman/message/30435186/

rule subset_vcf:
    input:
        vcfin = "data/processed/phased/{chr}.phased.filtered.vcf.gz"
    output:
        vcfout = "models/msmc/vcf/{sample}.{chr}.phased.vcf.gz"
    run:
        shell("bcftools view -s {wildcards.sample} -o {output.vcfout} -O z {input.vcfin}")

rule get_depth:
    input:
        bam = "data/raw/sorted_reads/{sample}.sorted.bam"
    output:
        depth = temp("models/msmc/indiv_masks/{sample}.{chr}.depth")
    run:
        shell("samtools depth -r {wildcards.chr} {input.bam} | \
        awk '{{sum += $3}} END {{print sum / NR}}' > {output.depth}")

rule indiv_mask:
    input:
        vcf = "models/msmc/vcf/{sample}.{chr}.phased.vcf.gz",
        depth = "models/msmc/indiv_masks/{sample}.{chr}.depth"
    output:
        mask = "models/msmc/indiv_masks/{sample}.{chr}.mask.bed.gz"
    params:
        mask = "models/msmc/indiv_masks/{sample}.{chr}.mask.bed"
    run:
        import subprocess
        depth = open(input.depth, "r")
        callString = "./src/MSMC_Mask.py -v " + str(input.vcf) + \
        " -o " + str(params.mask) + \
        " -d " + depth.read().rstrip() + " -g"
        print(callString)
        shell("python3 " + callString)
        shell("bgzip {params.mask}")


# generate input files for msmc
### NEED TO EXPAND INPUTS TO LIST ALL SAMPLES WITH FLAGS
rule msmc_input:
    input:
        expand("models/msmc/indiv_masks/{sample}.{{chr}}.mask.bed.gz", sample = SAMP_MSMC),
    output:
        "models/msmc/input/capitata.{chr}.multihetsep.txt"
    params:
        map_mask = "data/processed/mappability_masks/Boleracea_chr{chr}.mask.bed.gz",
        vcf = expand("models/msmc/vcf/{sample}.{{chr}}.phased.vcf.gz", sample = SAMP_MSMC),
        masks = lambda wildcards, input: " --mask=".join(input)
    shell:
        "./{msmc_tools}/generate_multihetsep.py \
        --mask={params.masks} \
        --mask={params.map_mask} \
        {params.vcf} > {output}"
