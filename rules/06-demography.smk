# TODO: edit input/output files - this is a positive mask that needs to be converted


# Demographic history with SMC++
# adapted from lovely script by cattlefriends Harly Durbin & Troy Rowan
# Check out SMC++ git repo for additional documentation:
#   https://github.com/popgenmethods/smcpp

################################################################################
# STEP 1: mask repetitive regions (e.g. centromeres)
#   SNPable program
#   generates a mask for single-end reads of length 'k' and stringency 'r'
#   default is k = 35
#   Thanks to Li Wang for helpful code!
#   Code also adapted from https://github.com/JackyHess/MSMC_analysis
#   splits reference into overlapping kmer sequences (100 bp)
#   the first step splits the reference fasta into read lengths of 100 bp
#   the pipe combines -l reads into a short read file
#   aligns reads to the genome using bwa aln
#   generates sam file for single end reads using bwa samse
#   perl script from SNPable to generate the raw (fasta) mask (gen_raw_mask.pl)
#   use gen_mask function from SNPable to apply stringency (-r)
#   use makeMappabilityMask.py to convert fasta to bed
#   requires download of makeMappabilityMask.py from https://github.com/stschiff/msmc-tools
# NOTE: had to update gzip.open in makeMappabilityMask.py from 'w' to 'wt'
################################################################################

rule make_mask:
    input:
        "data/external/ref/Boleracea_chromosomes.fasta"
    output:
        "data/processed/masks/Boleracea_{chr}.mask.bed.gz"
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


################################################################################
# STEP 2
# Create a list of distinguished individuals to use for SMC++ input files
# These are 5 randomly sampled individuals from each population
# Can be used to estimate a composite likelihood and uncertainty in the recent past
################################################################################

# rule distinguished_inds:
#     input:
#         pop_ids = config['pop_ids']
#     output:
#         "models/smc/distinguished_individuals.txt"
#     run:
#         import pandas as pd
#         pop_ids = pd.read_table("models/sample_ids.txt", names = ['Sample_ID', 'Population'])
#         smc_df = pop_ids.groupby(['Population']).filter(lambda x: len(x) > 4)
#         distind_ids = smc_df.groupby('Population', as_index = False).apply(pd.DataFrame.sample, n = 5, replace = False)
#         distind_ids.to_csv('models/smc/distinguished_individuals.txt', header = True, index = False, sep = '\t')

################################################################################
# STEP 3
# Create input files for SMC++ (.smc.gz format)
#   --mask = poorly mapped regions to exclude from analysis (otherwise assumes missing data = regions of homozygosity)
################################################################################

rule vcf2smc:
    input:
        vcf = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf",
        index = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.idx",
        mask = "data/processed/mappability_masks/scratch/{chr}.mask.bed.gz",
        dist_ind = "models/smc/distinguished_individuals.txt"
    output:
        "models/smc/input/{population}.{distinguished_ind}.{chr}.smc.gz"
    params:
        chrom = "{chr}",
    	# pop_choose exports a string of all individuals in current population
        list = pop_choose,
        # specify a distinguished individual
        # iterates over 10 randomly sampled individuals per population
        distind = "{distinguished_ind}"
    singularity:
        "docker://terhorst/smcpp:latest"
    run:
    # argument order:
    # mask, input vcf, output directory, chromosome, population: list of individuals
        shell("bgzip {input.vcf}")
        shell("tabix -p vcf {input.vcf}.gz")
        shell("smc++ vcf2smc \
        --mask {input.mask} \
        -d {params.distind} {params.distind} \
        {input.vcf}.gz {output} {params.chrom} {params.list}")
