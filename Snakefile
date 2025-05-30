################################################################################
## SNP calling pipeline for Brassica oleracea
## Follows GATK best practices
## https://software.broadinstitute.org/gatk/documentation/article.php?id=3893
## Sarah Turner-Hissong
################################################################################

configfile: "config.yaml"
include: "rules/global.smk"

# requires a conda environment to run:
# 	conda env create --name bo-demography --file environment.yaml

# submit.sh - bash script to run workflow
# submit.json - configuration file for cluster settings

################################################################################
##  a pseudo-rule that collects the target files (expected outputs)
################################################################################

rule all:
	input:
	# # SEQUENCE QUALITY
	# 	# fastqc = expand("qc/fastqc/{sample}_{readgroup}_fastqc.zip", sample = SAMPLES2, readgroup = ['R1', 'R2']),
	# 	multiqc = expand("qc/STJRI0{lane}_multiqc.html", lane = [1,2,3]),
	# # MAPPING
	# 	get_ref = "data/external/ref/Boleracea_chromosomes.fasta",
	# 	sort_bam = expand("data/raw/sorted_reads/{sample}.sorted.bam", sample = SAMPLES),
	# 	bamqc = expand("reports/bamqc/{sample}_stats/qualimapReport.html", sample = SAMPLES),
	# 	multibamqc = "reports/multisampleBamQcReport.html",
	# # CALLING
	# 	hap_caller = expand("data/interim/gvcf_files_bpres/{sample}.raw.snps.indels.g.vcf", sample = SAMPLES),
	# 	split_intervals = expand("data/processed/scattered_intervals/{intervals}-scattered.intervals",
	# 	intervals = INTERVALS),
	# 	joint_geno = expand("data/raw/vcf_bpres/{intervals}.raw.snps.indels.vcf", intervals = INTERVALS),
	# # FILTERING
	# 	filtered_invariant = expand("data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.invariant.sites.vcf", chr = CHR),
	# 	merge_filtered_snps = "data/processed/filtered_vcf_bpres/allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz",
	# 	# merge_filtered_allsites = expand("data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.allsites.vcf.gz", chr = CHR),
	# 	ind_stats = "reports/filtered.qual.dp5_200.maxnocall10.imiss",
	# # POPULATION STRUCTURE
	# 	admix_input = "data/processed/biallelic_snps.geno10.maf05.ldpruned.bed",
	# 	admixture = expand("models/admixture/biallelic_snps.geno10.maf05.ldpruned.{k}.Q", k = range(2, 15)),
	# 	pixy_pi = expand("models/pixy/B_oleracea_grouped_{chr}_{window_size}bp_{stat}.txt", chr = CHR, window_size = [10000, 50000, 100000], stat = ['pi', 'dxy', 'fst']),
	# SELECTIVE SWEEPS
		# create_bed = expand("models/RAiSD/{chr}_excluded_regions.bed", chr = CHR),
		raisd = expand("models/RAiSD/RAiSD_Report.{population}_{window_size}bp.{chr}", population = pop_dict.keys(), chr = CHR, window_size = [20,50]),
	# # DEMOGRAPHY
	# 	vcf2smc = smc_input_files,
	# 	smc_cv = expand("models/smc_cv_no_timepoints/{population}/model.final.json", population = distind_dict.keys()),
	# 	plot_estimate = "reports/smc_cv_no_timepoints_results.png",
	# 	smc_bootstrap = smc_bootstrap_input,
	# 	smc_cv_bootstrap = expand("models/smc_cv_bootstrap/{population}_{n_bootstrap}/model.final.json", population = distind_dict.keys(), n_bootstrap = range(1,11)),
	# 	joint_vcf2smc12 = smc_split_input_files12,
    #     joint_vcf2smc21 = smc_split_input_files21,
	# 	smc_split = expand("models/smc_split/{pop_pair}/model.final.json", pop_pair = pop_pair_dict.keys()),
		# plot_split = "reports/smc_split/brassica_all_pops_smc_split.png",
	# DADI
		# build_sfs_subpops = expand("models/dadi/sfs/{model}.fs", model = ['cap_gem_vir', 'sab_palm_alb', 'ital_botr']),
		# build_sfs_wild_dom = expand("models/dadi/sfs/{model}.fs", model = ['wild_domesticated']),
		# build_sfs_kales = expand("models/dadi/sfs/{model}.fs", model = ['wild_kale', 'gong_ital_kale']),
		# run_dadi_inference = expand("models/dadi/results/{model}.csv", model = ['cap_gem_vir', 'sab_palm_alb', 'ital_botr', 'wild_domesticated', 'wild_kale', 'gong_ital_kale'])


################################################################################
## Rule files to include
################################################################################

# include: "rules/01-fastqc.smk"
# include: "rules/02-mapping.smk"
# include: "rules/03-calling.smk"
# include: "rules/04-filtering.smk"
include: "rules/05-pop_structure.smk"
include: "rules/06-selective_sweeps.smk"
# include: "rules/07-demography.smk"
# include: "rules/08-split-times.smk"
include: "rules/09-dadi.smk"

