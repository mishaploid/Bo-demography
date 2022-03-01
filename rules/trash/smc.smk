# Demographic history with SMC++
# adapted from lovely script by cattlefriends Harly Durbin & Troy Rowan

# This command will parse data for the contig chr1 for samples S1 and S2
# which are members of population Pop1
# produces output for each chromosome

rule marginal_vcf2smc:
    input:
        vcf = "data/processed/filtered_snps/{chr}.filtered.dp6_200.nocall.snps.vcf.gz",
        index = "data/processed/filtered_snps/{chr}.filtered.dp6_200.nocall.snps.vcf.gz.tbi"
        # mask = "data/processed/masks/Boleracea_chr{chr}.mask.bed.gz"
    params:
        chrom = "{chr}",
		#For the current population, string of all individuals in that population
        list = pop_choose,
        mask = "data/processed/mappability_masks/scratch/{chr}.mask.bed.gz"
    output:
        "models/smc/input/{pop}.{chr}.smc.gz"
    shell:
	#Specify input vcf, specify output directory, speficy chromosome, specify string of list of individuals in the population
        "smc++ vcf2smc --mask {params.mask} {input.vcf} {output} {params.chrom} {params.list}"

# rule marginal_estimate:
# 	input:
# 		smc_chr = lambda wildcards: expand('models/smc/input/{pop}.{chr}.smc.gz', chr = chr, pop = wildcards.pop)
# 	output:
# 		"models/smc/estimate/{pop}/model.final.json"
# 	benchmark: "benchmarks/estimate/{pop}.estimate.benchmark"
# 	threads: 18
# 	params:
# 		model_in = "models/smc/input/{pop}.*",
# 		model_out_dir = "models/smc/estimate/{pop}",
# 		mu = 7e-9
# 	run:
# 		shell("smc++ estimate --cores 18 -o {params.model_out_dir} {params.mu} {params.model_in}")
# 		shell("smc++ estimate --cores 18 --timepoints 1e2 1e6 -o {params.model_out_dir}.timepoints {params.mu} {params.model_in}")

rule smc_cv:
	input:
		smc_chr = lambda wildcards: expand('models/smc/input/{pop}.{chr}.smc.gz', chr = chr, pop = wildcards.pop)
	output:
		expand("models/smc/cv/{{pop}}/fold{fold}/model.final.json", fold = ['0','1'])
	threads: 16
	params:
		model_in = "models/smc/input/{pop}.*",
		model_out_dir = "models/smc/cv/{pop}",
		mu = 7e-9
	run:
		shell("smc++ cv \
        --cores {threads} \
        --spline cubic \
		-o {params.model_out_dir} {params.mu} {params.model_in}")
        #         --timepoints 1e3 1e6 \

rule smc_estimate:
	input:
		smc_chr = lambda wildcards: expand('models/smc/input/{pop}.{chr}.smc.gz', chr = chr, pop = wildcards.pop)
	output:
		"models/smc/estimate/{pop}/model.final.json"
	threads: 24
	params:
		model_in = "models/smc/input/{pop}.*",
		model_out_dir = "models/smc/estimate/{pop}",
		mu = 7e-9
	run:
		shell("smc++ estimate \
        --cores {threads} \
        --spline cubic \
		-o {params.model_out_dir} {params.mu} {params.model_in}")
#         --timepoints 1e2 1e6 \

# # #Generate vcf2smc files containing the joint frequency spectrum for both populations
# rule joint_vcf2smc:
#     input:
#         vcf = "data/processed/filtered_snps/{chr}.filtered.snps.vcf.gz",
#         index = "data/processed/filtered_snps/{chr}.filtered.snps.vcf.gz.tbi"
#     output:
#         pop_pair_out = "models/smc/split/{pop_pair}.{chr}.smc.gz"
#     threads: 12
#     params:
#     	chrom = "{chr}",
#         mask = "data/processed/mappability_masks/scratch/{chr}.mask.bed.gz",
#     	pop_pair_string = pair_string_choose
#     shell:
#     	"smc++ vcf2smc \
#         --cores {threads} \
#         -m {params.mask} \
#         {input.vcf} \
#         {output.pop_pair_out} {params.chrom} {params.pop_pair_string}"
#
# rule split:
#     input:
#         expand("models/smc/split/{pop_pair}.{chr}.smc.gz", pop_pair = ['botrytis_italica', 'italica_botrytis'], chr = chr),
#         expand("models/smc/input/{pop}.{chr}.smc.gz", pop = pops, chr = chr)
#     threads: 16
#     params:
#         model_out_dir = "models/smc/split/test",
#         marginal_models = model_chooser
#     output:
#         model_out = "models/smc/split/model.final.json"
#     shell:
#         "smc++ split \
#         -o {params.model_out_dir} \
#         --cores {threads} \
#         {params.marginal_models} \
#         {input}"


# rule smc_cv_time:
# 	input:
# 		smc_chr = lambda wildcards: expand('models/smc/input_c150k/{pop}.{chr}.smc.gz', chr = chr, pop = wildcards.pop)
# 	output:
# 		expand("models/smc/cv_time_c150k/{{pop}}/fold{fold}/model.final.json", fold = ['0','1'])
# 	threads: 16
# 	params:
# 		model_in = "models/smc/input_c150k/{pop}.*",
# 		model_out_dir = "models/smc/cv_time_c150k/{pop}",
# 		mu = 7e-9
# 	run:
# 		shell("smc++ cv --cores {threads} \
#         --timepoints 1e3 2e6 \
# 		-o {params.model_out_dir} {params.mu} {params.model_in}")
#
# # rule plot_estimate:
# # 	input:
# # 		model = expand('models/smc/estimate/{pop}/model.final.json', pop = pops)
# # 	output:
# # 		plot = "reports/figures/smc/missing.png"
# # 	shell:
# # 		"smc++ plot -c {output.plot} {input.model}"
#
#
# # rule plot_split:
# #     input:
# #         model = "models/split/{rundate}.{dataset}/{pop_pair}/model.final.json"
# #     output:
# #         plot = "reports/figures/{rundate}.{dataset}.{pop_pair}.split.png"
# #     shell:
# #         "smc++ plot -c -g 2 {output.plot} {input.model}"
