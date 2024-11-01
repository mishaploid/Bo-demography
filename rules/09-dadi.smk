# Demographic modeling with dadi

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
	shell:
		"python3 src/dadi/build_frequency_spectra.py \
    	--vcf {input.vcf} \
    	--pop_info {input.pop_file} \
		--model_config_file {input.model_config} \
		--model {params.model} \
    	--subsample"

################################################################################
# STEP 2: run demographic models
################################################################################

rule run_dadi_inference:
	input:
		sfs = "models/dadi/sfs/{model}.fs"
	output:
		dadi_results = "models/dadi/results/{model}.csv"
	params:
		model = "{model}"
	shell:
		"python3 src/dadi/run_inference.py \
		--model {params.model}"