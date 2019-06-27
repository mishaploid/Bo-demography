# Demographic history with SMC++
# adapted from lovely script by cattlefriends Harly Durbin & Troy Rowan

# This command will parse data for the contig chr1 for samples S1 and S2
# which are members of population Pop1
# produces output for each chromosome

rule marginal_vcf2smc:
    input:
        vcf = "data/processed/{chr}.filtered.snps.vcf.gz",
        index = "data/processed/{chr}.filtered.snps.vcf.gz.tbi",
        mask = "data/processed/masks/Boleracea_chr{chr}.mask.bed.gz"
    params:
        chrom = "{chr}",
		#For the current population, string of all individuals in that population
        list = pop_choose
    output:
        "models/smc/input/{rundate}.{dataset}.{pop}.{chr}.smc.gz"
    shell:
	#Specify input vcf, specify output directory, speficy chromosome, specify string of list of individuals in the population
        "smc++ vcf2smc -m {input.mask} {input.vcf} {output} {params.chrom} {params.list}"

rule marginal_estimate:
    input:
	#Need to use lambda function in order to correctly iterate over individuals in individualdict
        smc_chr = lambda wildcards: expand('models/smc/input/{rundate}.{dataset}.{pop}.{chr}.smc.gz',
        rundate = rundate,
		chr = chr,
        dataset = dataset,
        pop = wildcards.pop)
    benchmark:
        "benchmarks/estimate/{rundate}.{dataset}.{pop}.estimate.benchmark"
    threads: 18
    params:
        model_in = "models/smc/input/{rundate}.{dataset}.{pop}.*",
        model_out_dir = "models/smc/cv/{rundate}.{dataset}.{pop}",
        mu = 7e-9
    output:
        "models/smc/cv/{rundate}.{dataset}.{pop}/model.final.json"
    shell:
        "smc++ cv --cores 18 -o {params.model_out_dir} {params.mu} {params.model_in}"

rule plot_estimate:
    input:
        model = lambda wildcards: expand('models/smc/cv/{rundate}.{dataset}.{pop}/model.final.json', rundate = rundate, dataset = dataset, pop = wildcards.pop)
    output:
        plot = "reports/figures/smc/{rundate}.{dataset}.{pop}.png"
    shell:
        "smc++ plot -c -g 2 {output.plot} {input.model}"

# #Generate vcf2smc files containing the joint frequency spectrum for both populations
# rule joint_vcf2smc:
#     input:
#         vcf = "data/raw/smc_input.vcf.gz",
#         index = "data/raw/smc_input.vcf.gz.tbi",
#         mask = "data/raw/mask_ref.bed.gz"
#     benchmark:
#         "benchmarks/joint_vcf2smc/{rundate}.{dataset}.{pop_pair}.{chr}.joint_vcf2smc.benchmark"
#     threads: 2
# 	params:
# 		chrom = "{chr}",
# 		pop_pair_string = pair_string_choose
# 	output:
# 		pop_pair_out = "data/processed/joint_vcf2smc/{rundate}.{dataset}.{pop_pair}.{chr}.smc.gz"
# 	shell:
# 		"smc++ vcf2smc --cores 8 -m {input.mask} {input.vcf} {output.pop_pair_out} {params.chrom} {params.pop_pair_string}"
#
# rule split:
#     input:
#         joint_vcf2smc = expand("data/processed/joint_vcf2smc/{{rundate}}.{{dataset}}.{{pop_pair}}.{chr}.smc.gz", chr = chr)
#     benchmark:
#         "benchmarks/split/{rundate}.{dataset}.{pop_pair}.split.benchmark"
#     threads: 8
#     params:
#         model_out_dir = "models/split/{rundate}.{dataset}/{pop_pair}",
#         marginal_models = model_chooser
#     output:
#         model_out = "models/split/{rundate}.{dataset}/{pop_pair}/model.final.json"
#     shell: "smc++ split -o {params.model_out_dir} --cores 8 {params.marginal_models} {input.joint_vcf2smc}"
#
# rule plot_split:
#     input:
#         model = "models/split/{rundate}.{dataset}/{pop_pair}/model.final.json"
#     output:
#         plot = "reports/figures/{rundate}.{dataset}.{pop_pair}.split.png"
#     shell:
#         "smc++ plot -c -g 2 {output.plot} {input.model}"
