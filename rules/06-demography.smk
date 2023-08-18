# TODO: edit input/output files - this is a positive mask that needs to be converted


# Demographic history with SMC++
# adapted from lovely script by cattlefriends Harly Durbin & Troy Rowan
# Check out SMC++ git repo for additional documentation:
#   https://github.com/popgenmethods/smcpp

################################################################################
# STEP 1: mask repetitive regions (e.g. centromeres) with SNPable program
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
#   NOTE: had to update gzip.open in makeMappabilityMask.py from 'w' to 'wt'
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
#         pop_ids = config['pruned_pop_ids']
#     output:
#         "models/smc/distinguished_individuals.txt"
#     run:
#         import pandas as pd
#         pop_ids = pd.read_table("models/pruned_sample_ids.txt", names = ['Sample_ID', 'Population'])
#         smc_df = pop_ids.groupby(['Population']).filter(lambda x: len(x) > 4)
#         distind_ids = smc_df.groupby('Population', as_index = False).apply(pd.DataFrame.sample, n = 5, replace = False)
#         distind_ids.to_csv('smc/distinguished_individuals.txt', header = True, index = False, sep = '\t')

################################################################################
# STEP 3
# Create input files for SMC++ (.smc.gz format)
#   --mask = poorly mapped regions to exclude from analysis (otherwise assumes missing data = regions of homozygosity)
################################################################################

rule vcf2smc:
    input:
        vcf = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz",
        index = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz.tbi",
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
    shell:
    # argument order:
    # mask, input vcf, output directory, chromosome, population: list of individuals
        "smc++ vcf2smc \
        --mask {input.mask} \
        -d {params.distind} {params.distind} \
        {input.vcf} {output} {params.chrom} {params.list}"


###############################################################################
# STEP 4
# Fit population size history to data
#   cv method uses cross-validation to obtain model parameters
################################################################################

rule smc_cv:
    input:
    	smc_cv_input
    output:
    	# cv_folds = expand("models/smc_cv_no_timepoints/{{population}}/fold{fold}/model.final.json", fold = ['0','1']),
        final_model = "models/smc_cv_no_timepoints/{population}/model.final.json"
    threads: 16
    params:
    	model_in = "models/smc/input/{population}.*",
    	model_out_dir = "models/smc_cv_no_timepoints/{population}",
    	mu = config['mu']
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
    	"smc++ cv \
        --cores {threads} \
        --spline cubic \
    	-o {params.model_out_dir} {params.mu} {params.model_in}"

################################################################################
# STEP 5
# Plot results for population size history estimates
#   --csv exports a csv formatted file with results
#   -g specifies the number of years per generation
################################################################################

rule plot_estimate:
    input:
        smc_out = expand("models/smc_cv_no_timepoints/{population}/model.final.json", population = distind_dict.keys())
    output:
        "reports/smc_cv_no_timepoints_results.png"
    params:
        gen = config['gen']
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
        "smc++ plot \
        --csv \
        -g {params.gen} \
        {output} \
        {input.smc_out}"

################################################################################
# STEP 6
# Bootstrap input files to get estimate of uncertainty in recent past
################################################################################

rule smc_bootstrap_input:
    input:
        expand("models/smc/input/{{population}}.{{distinguished_ind}}.{chr}.smc.gz", chr = CHR)
    output:
        expand("models/smc/bootstrap_input/{{population}}.{{distinguished_ind}}_rep_{n_bootstrap}/bootstrap_chr{boot_chr}.gz", n_bootstrap = range(1,11), boot_chr = range(1,10))
    params:
        population = "{population}",
        distind = "{distinguished_ind}",
        nr_bootstraps = 10,
        chunk_size = 5000000,
        nr_chr = 9,
        input_dir = "models/smc/input/{population}.{distinguished_ind}*"
    shell:
        "python3 src/smc_bootstrap.py \
        --nr_bootstraps {params.nr_bootstraps} \
        --chunk_size {params.chunk_size} \
        --chunks_per_chromosome 10 \
        --nr_chromosomes {params.nr_chr} \
        models/smc/bootstrap_input/{params.population}.{params.distind}_rep \
        {params.input_dir}"

rule smc_cv_bootstrap:
    input:
        smc_cv_boot_input
    output:
    	cv_folds = expand("models/smc_cv_bootstrap/{{population}}_{{n_bootstrap}}/fold{fold}/model.final.json", fold = ['0','1']),
        final_model = "models/smc_cv_bootstrap/{population}_{n_bootstrap}/model.final.json"
    threads: 16
    params:
    	model_in = "models/smc/bootstrap_input/{population}*rep_{n_bootstrap}/*",
    	model_out_dir = "models/smc_cv_bootstrap/{population}_{n_bootstrap}",
    	mu = config['mu']
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
    	"smc++ cv \
        --cores {threads} \
        --spline cubic \
    	-o {params.model_out_dir} {params.mu} {params.model_in}"

rule plot_bootstrap:
    input:
        boot_out = expand("models/smc_cv_bootstrap/{population}_{n_bootstrap}/model.final.json", population = distind_dict.keys(), n_bootstrap = range(1,11))
    output:
        "reports/smc_cv_no_timepoints_bootstrap_results.png"
    params:
        gen = config['gen']
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
        "smc++ plot \
        --csv \
        -g {params.gen} \
        {output} \
        {input.smc_out}"
<<<<<<< HEAD

# ################################################################################
# # STEP 7
# # Create input files for smc++ split command
# # joint site frequency spectrum for two populations
# ################################################################################
#
# rule joint_vcf2smc12:
#     input:
#         vcf = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz",
#         index = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz.tbi",
#         mask = "data/processed/mappability_masks/scratch/{chr}.mask.bed.gz"
#     output:
#         out12 = "models/smc_split_input/{pop_pair}_12.{distinguished_ind1}.{chr}.smc.gz",
#     params:
#     	chrom = "{chr}",
#         distind1 = "{distinguished_ind1}",
#     	pop_pair_string12 = pair_string_choose12
#     singularity:
#         "docker://terhorst/smcpp:latest"
#     shell:
#         """
#         smc++ vcf2smc \
#         --mask {input.mask} \
#         -d {params.distind1} {params.distind1} \
#         {input.vcf} {output.out12} {params.chrom} {params.pop_pair_string12}
#         """
#
# rule joint_vcf2smc21:
#     input:
#         vcf = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz",
#         index = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz.tbi",
#         mask = "data/processed/mappability_masks/scratch/{chr}.mask.bed.gz"
#     output:
#         out21 = "models/smc_split_input/{pop_pair}_21.{distinguished_ind2}.{chr}.smc.gz"
#     params:
#     	chrom = "{chr}",
#         distind2 = "{distinguished_ind2}",
#         pop_pair_string21 = pair_string_choose21
#     singularity:
#         "docker://terhorst/smcpp:latest"
#     shell:
#         """
#         smc++ vcf2smc \
#         --mask {input.mask} \
#         -d {params.distind2} {params.distind2} \
#         {input.vcf} {output.out21} {params.chrom} {params.pop_pair_string21}
#         """
=======
>>>>>>> 86e618b2156c55942bad0b79fb75df8e578cc3a0
