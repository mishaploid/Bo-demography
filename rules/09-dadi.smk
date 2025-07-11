# Demographic modeling with dadi

# Times in model are not estimating time from the present that something happened 
# estimating the length of a time interval 

################################################################################
# STEP 1: generate joint SFS for population comparisons
################################################################################

rule build_sfs_subpops:
	input: 
		vcf = "data/processed/filtered_vcf_bpres/allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz",
		pop_file = "models/pruned_sample_ids.txt",
		model_config = "src/dadi/model_config_subpops.json"
	output:
		joint_sfs = "models/dadi/sfs/{model}.fs"
	params:
		model = "{model}"
	wildcard_constraints:
		model = "cap_gem_vir|sab_palm_alb|ital_botr"
	resources:
		mem_mb = 50000
	conda:
		"bo-demography"
	shell:
		"python3 src/dadi/build_frequency_spectra.py \
    	--vcf {input.vcf} \
    	--pop_info {input.pop_file} \
		--model_config_file {input.model_config} \
		--model {params.model} \
    	--subsample"

rule build_sfs_wild_dom:
	input: 
		vcf = "data/processed/filtered_vcf_bpres/allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz",
		pop_file = "models/pruned_sample_ids_wild_dom.txt",
		model_config = "src/dadi/model_config_subpops.json"
	output:
		joint_sfs = "models/dadi/sfs/{model}.fs"
	params:
		model = "{model}"
	wildcard_constraints:
		model = "wild_domesticated"
	resources:
		mem_mb = 50000
	conda:
		"../environment.yaml"
	shell:
		"python3 src/dadi/build_frequency_spectra.py \
    	--vcf {input.vcf} \
    	--pop_info {input.pop_file} \
		--model_config_file {input.model_config} \
		--model {params.model} \
    	--subsample"

rule build_sfs_kales:
	input: 
		vcf = "data/processed/filtered_vcf_bpres/allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz",
		pop_file = "models/pruned_sample_ids_wild_kale.txt",
		model_config = "src/dadi/model_config_subpops.json"
	output:
		joint_sfs = "models/dadi/sfs/{model}.fs"
	params:
		model = "{model}"
	wildcard_constraints:
		model = "wild_kale|gong_ital_kale"
	resources:
		mem_mb = 50000
	conda:
		"../environment.yaml"
	shell:
		"python3 src/dadi/build_frequency_spectra.py \
    	--vcf {input.vcf} \
    	--pop_info {input.pop_file} \
		--model_config_file {input.model_config} \
		--model {params.model} \
    	--subsample"

################################################################################
# STEP 2: run demographic model optimization 
################################################################################

rule run_dadi_inference:
	input:
		sfs = "models/dadi/sfs/{model}.fs",
		model_config = "src/dadi/model_config_subpops.json"
	output:
		dadi_results = "models/dadi/results/{model}.csv"
	params:
		model = "{model}"
	resources:
		mem_mb = 50000
	conda:
		"../environment.yaml"
	shell:
		"python3 src/dadi/run_inference.py \
		--model_config_file {input.model_config} \
		--model {params.model}"


# read in csv file, sort by likelihood and select best option (least negative number) 
# values are in genetic units, need to convert to real units using length of sequence 
# assumed mutation rate, etc. 
# sequence length can be size of genome (should remove masked regions)
# uncertainty estimates around parameters (built in methods in dadi) 
# for bootstrapping: dadi function that will take VCF, resample, and create an SFS 
#	do not need to rerun inference, will use ML estimate of parameters to look at uncertainty and likelihood around bootstrapped SFS 
#	if run time is prohibitive, can try a different function to create data dictionary from VCF 
# def bootstraps_subsample_vcf
# (
# vcf_filename, popinfo_filename, subsample, Nboot, chunk_size, pop_ids, filter=True, flanking_info=[None, None], mask_corners=True, polarized=True)
# maximum likelihood parameter estimates + boostrapped samples of SFS 
# uses ML parameter estimates and bootstrapped SFS to explore likelihood surface around best estimates to acquire confidence intervals 
# approximates properties of likelihood surface for CIs without redoing inference 
