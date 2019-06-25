# Demographic history with SMC++
# adapted from lovely script by Harly Durbin & Troy Rowan

#Dictionary of sample id string for each per population
capitata = ['SamC_' + str(x).rjust(3, '0') for x in range(1,46)]
capitatastr = ",".join(capitata)
gongylodes = ['SamC_' + str(x).rjust(3, '0') for x in range(46,65)]
gongylodesstr = ",".join(gongylodes)
italica = ['SamC_' + str(x).rjust(3, '0') for x in range(85,108)]
italicastr = ",".join(italica)
botrytis = ['SamC_' + str(x).rjust(3, '0') for x in range(65,85)]
botrytisstr = ",".join(botrytis)
alboglabra = ['SamC_108', 'SamC_109', 'SamC_110', 'SamC_111']
alboglabrastr = ",".join(alboglabra)
gemmifera = ['SamC_112', 'SamC_113']
gemmiferastr = ",".join(gemmifera)
acephala = ['SamC_114', 'SamC_115']
acephalastr = ",".join(acephala)
sabellica = ['SamC_116', 'SamC_117']
sabellicastr = ",".join(sabellica)
wild = ['SamC_118', 'SamC_119']
wildstr = ",".join(wild)
kale = ['SamC_114', 'SamC_115', 'SamC_116', 'SamC_117', 'SamC_108', 'SamC_109', 'SamC_110', 'SamC_111']
kalestr = ",".join(kale)

popdict = {
'capitata':'capitata:' + capitatastr,
'gongylodes':'gongylodes:' + gongylodesstr,
'italica':'italica:' + italicastr,
'botrytis':'botrytis:' + botrytisstr,
'alboglabra' : 'alboglabra:' + alboglabrastr,
'gemmifera' : 'gemmifera:' + gemmiferastr,
'acephala' : 'acephala:' + acephalastr,
'sabellica' : 'sabellica:' + sabellicastr,
'wild' : 'wild:' + wildstr,
'kale' : 'kale:' + kalestr}

# Dictionary of population pairs for split
pop_pair_dict = {'botrytis_italica':['botrytis', 'italica'],
'italica_botrytis':['italica', 'botrytis'],
'italica_wild':['italica', 'wild'],
'wild_italica':['wild', 'italica']}

# --thinning controls frequency that CSFS is emitted
# decreasing can potentially increase accuracy of results in recent past
# Default is 1000*log(n), change it to 100*log(n) to provide better resolution in more recent past
# n = haploid sample size in undistinguished portion
# use 100*log((N*2)-1)
thinning_dict = {
"capitata": ['195'],
"gongylodes": ['157'],
"italica": ['165'],
"botrytis": ['159'],
"alboglabra": ['117'],
"gemmifera": ['60'],
"acephala": ['60'],
"sabellica": ['60'],
"wild": ['60'],
"kale": ['120']}

def pop_choose(WC):
	list = popdict[WC.pop]
	return list

def thinning_choose(WC):
	thinning_k = thinning_dict[WC.pop]
	return thinning_k

#Based on current wildcard, find pop pair in pop_pair_dict
#For each population in current pop pair, find pop in popdict
#Combine strings to create a super long string
#Return string
def pair_string_choose(WC):
	pops = pop_pair_dict[WC.pop_pair]
	pop1 = popdict[pops[0]]
	pop2 = popdict[pops[1]]
	pop_pair_string = pop1 + " " + pop2
	return pop_pair_string

#Based on current wildcard, find pop pair in pop_pair_dict
#Create string: location of model.final.json for each pop
def model_chooser(WC):
	pops = pop_pair_dict[WC.pop_pair]
	pop1 = pops[0]
	pop2 = pops[1]
	model_string = "models/estimate/" + WC.rundate + "." + WC.dataset + "/" + pop1 + "/model.final.json " + "models/estimate/" + WC.rundate + "." + WC.dataset + "/" + pop2 + "/model.final.json"
	return model_string

#"This command will parse data for the contig chr1 for samples S1 and S2 which are members of population Pop1. You should run this once for each independent contig in your dataset, producing one SMC++ output file per contig."
rule marginal_vcf2smc:
    input:
        vcf = "data/raw/smc_input.vcf.gz",
        index = "data/raw/smc_input.vcf.gz.tbi",
        mask = "data/raw/mask_ref.bed.gz"
    benchmark:
        "benchmarks/vcf2smc/190107.{dataset}.{pop}.{chr}.vcf2smc.benchmark"
    threads: 1
    params:
        chrom = "{chr}",
		#For the current population, string of all individuals in that population
        list = pop_choose
		#Select distinguished individual. Need to iterate over every chromosome of every possible distinguished individual
		# distinguished = "{distinguished}"
    output:
        smc = "data/processed/190107.{dataset}.{pop}.{chr}.smc.gz"
    shell:
	#Specify input vcf, specify output directory, speficy chromosome, specify string of list of individuals in the population
        "smc++ vcf2smc -m {input.mask} {input.vcf} {output.smc} {params.chrom} {params.list}"

# Generate population-specific estimates
# rule marginal_estimate:
#     input:
# 	#Need to use lambda function in order to correctly iterate over individuals in individualdict
#         smc_chr = lambda wildcards: expand('data/processed/190107.{dataset}.{pop}.{chr}.smc.gz',
#         chr = CHR,
#         dataset = wildcards.dataset,
#         pop = wildcards.pop)
#     benchmark:
#         "benchmarks/estimate/{rundate}.{dataset}.{pop}.estimate.benchmark"
#     threads: 8
#     params:
#         model_in = "data/processed/190107.{dataset}.{pop}.*",
#         model_out_dir = "estimate/{rundate}.{dataset}/{pop}",
#         mu = 7e-9
#     output:
#         model_plot = "estimate/{rundate}.{dataset}/{pop}/model.final.json",
#         model_iter = expand('estimate/{{rundate}}.{{dataset}}/{{pop}}/fold{fold}/.model.iter1.json', fold = [0,1])
# 	#specify -v for verbose to increase debuggling output
# 	#update: as of v1.5, estimate command renamed cv
# 	#use --spline cubic to get smoothed representation
#     shell:
#         "smc++ cv --cores 8 -o {params.model_out_dir} {params.mu} {params.model_in}"

rule timepoint_estimate:
    input:
	#Need to use lambda function in order to correctly iterate over individuals in individualdict
        smc_chr = lambda wildcards: expand('data/processed/190107.{dataset}.{pop}.{chr}.smc.gz', chr = CHR, dataset = wildcards.dataset, pop = wildcards.pop)
    benchmark:
        "benchmarks/estimate/{rundate}.{dataset}.{pop}.timepoints.estimate.benchmark"
    threads: 8
    params:
        model_in = "data/processed/190107.{dataset}.{pop}.*",
        model_out_dir = "models/estimate/{rundate}.{dataset}.timepoints/{pop}",
        mu = 7e-9,
        thinning_k = thinning_choose
    output:
        model = "models/estimate/{rundate}.{dataset}.timepoints/{pop}/model.final.json"
    shell:
        "smc++ estimate --cores 8 --timepoints 50 150000 --thinning {params.thinning_k} -o {params.model_out_dir} {params.mu} {input.smc_chr}"

### To run cv, change output and shell command
        # model_iter = expand('estimate/{{rundate}}.{{dataset}}.timepoints/{{pop}}/fold{fold}/.model.iter1.json', fold = [0,1])
        # "smc++ cv --cores 36 --timepoints 50 150000 --thinning {params.thinning_k} -o {params.model_out_dir} {params.mu} {params.model_in}"

rule plot_estimate:
    input:
        model = lambda wildcards: expand('models/estimate/{rundate}.{dataset}.timepoints/{pop}/model.final.json', rundate = wildcards.rundate, dataset = wildcards.dataset, pop = wildcards.pop)
    output:
        plot = "reports/figures/{rundate}.{dataset}.{pop}.timepoints.png"
    shell:
        "smc++ plot -c -g 2 {output.plot} {input.model}"

#Generate vcf2smc files containing the joint frequency spectrum for both populations
rule joint_vcf2smc:
    input:
        vcf = "data/raw/smc_input.vcf.gz",
        index = "data/raw/smc_input.vcf.gz.tbi",
        mask = "data/raw/mask_ref.bed.gz"
    benchmark:
        "benchmarks/joint_vcf2smc/{rundate}.{dataset}.{pop_pair}.{chr}.joint_vcf2smc.benchmark"
    threads: 2
	params:
		chrom = "{chr}",
		pop_pair_string = pair_string_choose
	output:
		pop_pair_out = "data/processed/joint_vcf2smc/{rundate}.{dataset}.{pop_pair}.{chr}.smc.gz"
	shell:
		"smc++ vcf2smc --cores 8 -m {input.mask} {input.vcf} {output.pop_pair_out} {params.chrom} {params.pop_pair_string}"

rule split:
    input:
        joint_vcf2smc = expand("data/processed/joint_vcf2smc/{{rundate}}.{{dataset}}.{{pop_pair}}.{chr}.smc.gz", chr = chr)
    benchmark:
        "benchmarks/split/{rundate}.{dataset}.{pop_pair}.split.benchmark"
    threads: 8
    params:
        model_out_dir = "models/split/{rundate}.{dataset}/{pop_pair}",
        marginal_models = model_chooser
    output:
        model_out = "models/split/{rundate}.{dataset}/{pop_pair}/model.final.json"
    shell: "smc++ split -o {params.model_out_dir} --cores 8 {params.marginal_models} {input.joint_vcf2smc}"

rule plot_split:
    input:
        model = "models/split/{rundate}.{dataset}/{pop_pair}/model.final.json"
    output:
        plot = "reports/figures/{rundate}.{dataset}.{pop_pair}.split.png"
    shell:
        "smc++ plot -c -g 2 {output.plot} {input.model}"
