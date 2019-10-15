################################################################################
## SNP calling pipeline for Brassica oleracea
## Follows GATK best practices
## https://software.broadinstitute.org/gatk/documentation/article.php?id=3893
## Sarah Turner-Hissong
################################################################################

import csv
import os

configfile: "config.yaml"

# requires a conda environment to run:
# conda env create --name bo-demography --file environment.yaml

# submit.sh - bash script to run workflow
# submit.json - configuration file for cluster settings

################################################################################
## Specify sample names
################################################################################

# Option 1: pull unique identifier from fastq files
# SAMPLES = glob_wildcards("data/external/fastq_raw/{sample}_pass_1.fastq.gz").sample

# Option 2: use names from sorted bams
# useful if you don't want to store fastq locally
SAMPLES = glob_wildcards("data/raw/sorted_reads/{sample}.sorted.bam").sample

# Option 3: manual entry of file names
# SAMPLES = ['SamC_' + str(x).rjust(3, '0') for x in range(1,120)]

# remove Brassica rapa for now
SAMPLES.remove('B_rapa')

################################################################################
## Dictionary of SRA identifiers for each sample
## useful if pulling entries from SRA
## download SRA project info first and edit filename in config.yaml
################################################################################

# create a dictionary of SRA identifiers from downloaded project info
with open(config['sra_info'], mode = 'r') as infile:
	reader = csv.reader(infile)
	sample_dict = {rows[29]:rows[0] for rows in reader}

# add Brassica rapa reference (Chiifu) short reads
sample_dict.update({'B_rapa':'SRR7881031'})

################################################################################
## List of chromosome names
################################################################################

chr = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']

################################################################################
##  a pseudo-rule that collects the target files (expected outputs)
################################################################################

rule all:
	input:
		# SEQUENCE QUALITY
		# fastqc = expand("qc/fastqc/{sample}_{readgroup}_fastqc.zip", sample = SAMPLES2, readgroup = ['R1', 'R2']),
		multiqc = expand("qc/STJRI0{lane}_multiqc.html", lane = [1,2,3]),
		# MAPPING
		get_ref = expand("data/external/ref/Boleracea_chromosomes.fasta"),
		sort_bam = expand("data/raw/sorted_reads/{sample}.sorted.bam", sample = SAMPLES),
		bamqc = expand("reports/bamqc/{sample}_stats/qualimapReport.html", sample = SAMPLES),
		multibamqc = "reports/multisampleBamQcReport.html",
		# CALLING
		hap_caller = expand("data/interim/gvcf_files/{sample}.raw.snps.indels.g.vcf", sample = SAMPLES),
		joint_geno = expand("data/raw/vcf/{chr}.raw.snps.indels.vcf", chr = chr),
		# FILTERING
		bgzip_vcf = expand("data/processed/filtered_snps/{chr}.filtered.snps.vcf.gz", chr = chr),
		merge_vcfs = "data/processed/filtered_snps/oleracea_combined.vcf.gz"


################################################################################
## Rule files to include
################################################################################

include: "rules/fastqc.smk"
include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/filtering.smk"
# include: "rules/msmc.smk"
# include: "rules/phasing.smk"
# include: "rules/smc.smk"
# include: "rules/demography.smk"
# include: "rules/admixture.smk"
# include: "rules/angsd.smk"

### intermediate steps (need to clean up & only include essential outputs in rule_all)

# admix_input
# "models/admixture/combined.pruned.bed",
# admixture
# expand("models/admixture/combined.pruned.{k}.Q", k = [1,2,3,4,5,7,8])
# # angsd_depth
# expand("reports/ALL.{chr}.qc.depthGlobal", chr = chr)
### SMC round 2
# smc = expand("models/smc/cv_1e3_1e6/{pop}/fold{fold}/model.final.json", pop = pops, fold = ['0','1'])
# smc = expand("models/smc/estimate_tp/{pop}/model.final.json", pop = pops),
# vcf2smc = expand("models/smc/input/{pop}.{chr}.smc.gz", pop = pops, chr = chr),
# jointvcf2smc = expand("models/smc/split/{pop_pair}.{chr}.smc.gz", pop_pair = ['wild_alboglabra', 'alboglabra_wild'], chr = chr)
# split = "models/smc/split/model.final.json"
### PHASING
# phased = expand("data/processed/phased/{chr}.phased.filtered.vcf.gz", chr = chr),
### MSMC
# depth = expand("models/msmc/indiv_masks/{sample}.{chr}.depth", sample = SAMP_MSMC, chr = chr),
# mappability = expand("data/processed/mappability_masks/Boleracea_chr{chr}.mask.bed.gz", chr = chr),
# masks = expand("models/msmc/indiv_masks/{sample}.{chr}.mask.bed", sample = SAMP_MSMC, chr = chr),
# phased = expand("models/msmc/vcf/{sample}.{chr}.phased.vcf.gz", sample = SAMP_MSMC, chr = chr),
# msmcin = expand("models/msmc/input/capitata.{chr}.multihetsep.txt", chr = chr)

################################################################################
# ## dictionaries for SMC++
# # create all of the dictionaries
# # sample ids for each population
# capitata = ','.join(['SamC_' + str(x).rjust(3, '0') for x in range(1,46)])
# gongylodes = ','.join(['SamC_' + str(x).rjust(3, '0') for x in range(46,65)])
# italica = ','.join(['SamC_' + str(x).rjust(3, '0') for x in range(85,108)])
# botrytis = ','.join(['SamC_' + str(x).rjust(3, '0') for x in range(65,85)])
# alboglabra = ','.join(['SamC_108', 'SamC_109', 'SamC_110', 'SamC_111'])
# gemmifera = ','.join(['SamC_112', 'SamC_113'])
# acephala = ','.join(['SamC_114', 'SamC_115'])
# sabellica = ','.join(['SamC_116', 'SamC_117'])
# wild = ','.join(['SamC_118', 'SamC_119'])
#
# pops = ['capitata', 'gongylodes', 'italica', 'botrytis', 'alboglabra', 'gemmifera', 'acephala', 'sabellica', 'wild']
#
# # pop_pair = ['botrytis_italica', 'italica_botrytis']
#
# # dictionary of population and corresponding samples
# popdict = {
# 'capitata':'capitata:' + capitata,
# 'gongylodes':'gongylodes:' + gongylodes,
# 'italica':'italica:' + italica,
# 'botrytis':'botrytis:' + botrytis,
# 'alboglabra' : 'alboglabra:' + alboglabra,
# 'gemmifera' : 'gemmifera:' + gemmifera,
# 'acephala' : 'acephala:' + acephala,
# 'sabellica' : 'sabellica:' + sabellica,
# 'wild' : 'wild:' + wild,
# }
#
# # dictionary of population pairs to estimate divergence
# pop_pair_dict = {'botrytis_italica':['botrytis', 'italica'],
# 'italica_botrytis':['italica', 'botrytis'],
# 'italica_wild':['italica', 'wild'],
# 'wild_italica':['wild', 'italica']}
#
# def pop_choose(WC):
# 	list = popdict[WC.pop]
# 	return list
#
# def pop_pair_choose(WC):
# 	list = popdict[WC.pop_pair]
# 	return list
#
# # for split time:
# # based on current wildcard, find pop pair in pop_pair_dict
# # for each population in current pop pair, find pop in popdict
# # combine strings to create a super long string
# # return the super long string
#
# def pair_string_choose(WC):
# 	pops = pop_pair_dict[WC.pop_pair]
# 	pop1 = popdict[pops[0]]
# 	pop2 = popdict[pops[1]]
# 	pop_pair_string = pop1 + " " + pop2
# 	return pop_pair_string
#
# # create string: location of model.final.json for each pop
#
# def model_chooser(WC):
# 	pops = pop_pair_dict[WC.pop_pair]
# 	pop1 = pops[0]
# 	pop2 = pops[1]
# 	model_string = "models/smc/cv_1e3_1e6/" + pop1 + "/model.final.json " \
# 	+ "models/smc/cv_1e3_1e6/" + pop2 + "/model.final.json"
# 	return model_string
