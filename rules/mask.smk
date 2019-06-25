# mask repetitive regions (e.g. centromeres)

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

rule make_mask:
    input:
        "data/external/ref/Boleracea_chromosomes.fasta"
    output:
        fa = "data/interim/masks/{chr}_mask_100_50.fa"
        # bed = expand("data/processed/masks/Boleracea_{chr}.mask.bed.gz", chr = chr)
    params:
        tmp = "/scratch/sdturner/{chr}",
        chr = "{chr}",
        newfasta = "data/external/ref/Boleracea_{chr}.fa"
    threads: 16
    run:
        shell("mkdir -p {params.tmp}")
        shell("samtools faidx {input} {params.chr} > {params.newfasta}")
        shell("{snpable}/splitfa {params.newfasta} 100 | split -l 80000000 - {params.tmp}/x")
        shell("cat {params.tmp}/x* > {params.tmp}/simu_reads.fq")
        shell("bwa aln -R 1000000 -O 3 -E 3 -t {threads} {input} {params.tmp}/simu_reads.fq > {params.tmp}/simu_reads.sai")
        shell("bwa samse {input} {params.tmp}/simu_reads.sai {params.tmp}/simu_reads.fq > {params.tmp}/simu_reads.sam")
        shell("perl {snpable}/gen_raw_mask.pl {params.tmp}/simu_reads.sam > {params.tmp}/raw_mask_100.fa")
        shell("{snpable}/gen_mask -l 100 -r 0.5  {params.tmp}/raw_mask_100.fa > {output.fa}")
        shell("rm -rf {params.tmp}")
        # shell("src/makeMappabilityMask.py {output.fa}")
