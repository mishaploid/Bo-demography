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
		model_config = "src/dadi/model_config.json"
	output:
		joint_sfs = "models/dadi/sfs/{model}.fs"
	params:
		model = "{model}"
	wildcard_constraints:
		model = "cap_gem_vir|sab_palm_alb|ital_botr"
	resources:
		mem_mb = 50000
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
		model_config = "src/dadi/model_config.json"
	output:
		joint_sfs = "models/dadi/sfs/{model}.fs"
	params:
		model = "{model}"
	wildcard_constraints:
		model = "wild_domesticated"
	resources:
		mem_mb = 50000
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
		model_config = "src/dadi/model_config.json"
	output:
		joint_sfs = "models/dadi/sfs/{model}.fs"
	params:
		model = "{model}"
	wildcard_constraints:
		model = "wild_kale|gong_ital_kale"
	resources:
		mem_mb = 50000
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

# added inbreeding coefficients as parameters
# plink2 --vcf ../../data/processed/filtered_vcf_bpres/allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz --het --out allsamps_inbreeding_coef --allow-extra-chr

rule run_dadi_inference:
	input:
		sfs = "models/dadi/sfs/{model}.fs",
		model_config = "src/dadi/model_config.json"
	output:
		dadi_results = "models/dadi/results/{model}.csv"
	params:
		model = "{model}"
	resources:
		mem_mb = 100000
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

################################################################################
# STEP 3: bootstrapping results to estimate uncertainty
################################################################################

# https://dadi.readthedocs.io/en/latest/api/dadi/Misc.html

# def bootstraps_subsample_vcf
# (
# vcf_filename, popinfo_filename, subsample, Nboot, chunk_size, pop_ids, filter=True, flanking_info=[None, None], mask_corners=True, polarized=True)

# inputs:
# 	vcf = the vcf (yay this is easy)
#	popinfo_filename = same file used for sfs 
# 	subsample = dictionary of number of individuals from each population
#	Nboot = number of bootstraps 
# 	chunk_size = use same size as for bootstrapping in SMC++ 
# 	pop_ids = which populations to sample for bootstrapping
#	filter = whether PASS/FAIL calls are removed from VCF 
#	flanking_info = 
#	mask_corners = whether to include or exclude monomorphic sites 
#	polarized = are calls using ancestral alleles as reference 

rule bootstrap_sfs_subpops:
	input:
		vcf = "data/processed/filtered_vcf_bpres/allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz",
		pop_file = "models/pruned_sample_ids.txt",
		model_config = "src/dadi/model_config.json"	
	output:
		bootstrapped_sfs = expand("models/dadi/sfs_bootstrap/{{model}}_{rep}.fs", rep = range(0,9))
	params:
		model = "{model}",
		rep = range(0,9)
	resources:
		mem_mb = 50000
	wildcard_constraints:
		model = "cap_gem_vir|sab_palm_alb|ital_botr"
	shell: 
		"python3 src/dadi/bootstrap_frequency_spectra.py \
    	--vcf {input.vcf} \
    	--pop_info {input.pop_file} \
		--model_config_file {input.model_config} \
		--model {params.model} \
		--subsample"

rule bootstrap_sfs_wild_domesticated:
	input:
		vcf = "data/processed/filtered_vcf_bpres/allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz",
		pop_file = "models/pruned_sample_ids_wild_dom.txt",
		model_config = "src/dadi/model_config.json"	
	output:
		bootstrapped_sfs = expand("models/dadi/sfs_bootstrap/{{model}}_{rep}.fs", rep = range(0,9))
	params:
		model = "{model}",
		rep = range(0,9)
	resources:
		mem_mb = 50000
	wildcard_constraints:
		model = "wild_domesticated"
	shell: 
		"python3 src/dadi/bootstrap_frequency_spectra.py \
    	--vcf {input.vcf} \
    	--pop_info {input.pop_file} \
		--model_config_file {input.model_config} \
		--model {params.model} \
		--subsample"

rule bootstrap_sfs_wild_kale:
	input:
		vcf = "data/processed/filtered_vcf_bpres/allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz",
		pop_file = "models/pruned_sample_ids_wild_kale.txt",
		model_config = "src/dadi/model_config.json"	
	output:
		bootstrapped_sfs = expand("models/dadi/sfs_bootstrap/{{model}}_{rep}.fs", rep = range(0,9))
	params:
		model = "{model}",
		rep = range(0,9)
	resources:
		mem_mb = 50000
	wildcard_constraints:
		model = "wild_kale|gong_ital_kale"
	shell: 
		"python3 src/dadi/bootstrap_frequency_spectra.py \
    	--vcf {input.vcf} \
    	--pop_info {input.pop_file} \
		--model_config_file {input.model_config} \
		--model {params.model} \
		--subsample"

# rule dadi_est_uncertainty:
	


# uncert = dadi.Godambe.GIM_uncert(func_ex, grid_pts, all_boot, p0, data, log, multinom, eps, return_GIM)
# func_ex = model config 
# grid_pts = same as previous
# all_boot = bootstrap output 
# Here func_ex is the model function, grid_pts is the set of grid points used in extrapolation, all_boot is a list containing bootstrapped data sets, p0 is the best-fit parameters, and data is the original data. If log = True, then uncertainties will be calculated for the logs of the parameters; these can be interpreted as relative uncertainties for the parameters themselves. If multinom = True, it is assumed that θ
#  is not an explicit parameter of the model (this is the most common case). eps is the relative step size to use when taking numerical derivatives; the default value is often sufficient. The returned uncert is an array equal in length to p0, where each entry in uncert is the estimated standard deviation of the parameter it corresponds to in p0. If multinom = True, there will be one extra entry in uncert, corresponding to θ
# . If return_GIM = True, then the return value will be (uncert, GIM), where GIM is the full Godambe Information Matrix, for use in propagating uncertainties.

