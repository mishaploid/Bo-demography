################################################################################
## SNP calling pipeline for Brassica oleracea
## Follows GATK best practices
## https://software.broadinstitute.org/gatk/documentation/article.php?id=3893
## Sarah Turner-Hissong
################################################################################

import csv
import os
import pandas as pd

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
SAMPLES_BAM = glob_wildcards("data/raw/sorted_reads/{sample}.sorted.bam").sample

# Option 3: manual entry of file names
# SAMPLES = ['SamC_' + str(x).rjust(3, '0') for x in range(1,120)]

# remove Brassica rapa for now
SAMPLES_BAM.remove('B_rapa')

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
# add additional Brassica cretica sequences from SRA
sample_dict.update({'B_rapa':'SRR7881031',
                    'B_cretica_A':'SRR9331103',
                    'B_cretica_B':'SRR9331104',
                    'B_cretica_C':'SRR9331105',
                    'B_cretica_D':'SRR9331106'})

# remove completed samples from SRA dictionary
SAMPLES_SRA = [*sample_dict] # unpack keys
SAMPLES_SRA.remove('SampleName')

for i in SAMPLES_BAM:
	try:
		SAMPLES_SRA.remove(i)
	except ValueError:
		pass

print(SAMPLES_SRA)

# combine sample names
SAMPLES = SAMPLES_BAM + SAMPLES_SRA

SAMPLES.remove('B_rapa')

################################################################################
## List of chromosome names
################################################################################

CHR = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']

################################################################################
## Number of intervals for GATK
## set to 200 intervals
################################################################################

INTERVALS = ["{:04d}".format(x) for x in list(range(200))]

################################################################################
## dictionaries for SMC++ cv command
# create all of the dictionaries
# sample ids for each population
# iteration over distinguished individuals
# https://stackoverflow.com/questions/18695605/python-pandas-dataframe-to-dictionary
################################################################################

# read in list of sample/morphotype ids
pop_ids = pd.read_csv(config['pop_ids'], sep="\t", names=['Sample_ID', 'Population'])

# filter for pops with more than 4 samples
pop_ids = pop_ids.groupby(['Population']).filter(lambda x: len(x) > 4)

# create a list of samples for each population
pop_list = pop_ids.groupby('Population')['Sample_ID'].apply(lambda x: x.tolist())

# convert sample list to a python dictionary object
pop_dict = pop_list.to_dict()

# join sample names and separate them with a ','
from collections import defaultdict
for key in pop_dict.keys():
    pop_dict[key] = ','.join(pop_dict[key])

# format input for SMC++
# Population:sample1,sample2,sample3...samplen
for key, value in pop_dict.items() :
    pop_dict[key] = key + ':' + value

# function to print formatted list of population: sample ids using wildcards
def pop_choose(wildcards):
	list = pop_dict[wildcards.population]
	return list

# Distinguished individuals to use for SMC++ input files
# These are 5 randomly sampled individuals from each population
# Can be used to estimate a composite likelihood and uncertainty in the recent past

distind_ids = pd.read_table(config['distinguished_individuals'])

distind_list = distind_ids.groupby('Population')['Sample_ID'].apply(lambda x: x.tolist())

distind_dict = distind_list.to_dict()

# create a list of filenames for smc++ input files using different distinguished individuals (output of vcf2smc)
smc_input_files = [expand('models/smc/input/{population}.{distinguished_ind}.{chr}.smc.gz',
                    population = key, distinguished_ind = value, chr = CHR)
                    for key, value in distind_dict.items()]

# create a list of filenames to feed into smc++ cv command
# want to include all .smc.gz files for each population (includes all chromosomes and distinguished individuals)
def smc_cv_input(wildcards):
    files = expand('models/smc/input/{population}.{distinguished_ind}.{chr}.smc.gz', chr=CHR, population=wildcards.population, distinguished_ind=distind_dict[wildcards.population])
    return files


################################################################################
##  a pseudo-rule that collects the target files (expected outputs)
################################################################################

rule all:
	input:
		# SEQUENCE QUALITY
		# fastqc = expand("qc/fastqc/{sample}_{readgroup}_fastqc.zip", sample = SAMPLES2, readgroup = ['R1', 'R2']),
		multiqc = expand("qc/STJRI0{lane}_multiqc.html", lane = [1,2,3]),
		# MAPPING
		get_ref = "data/external/ref/Boleracea_chromosomes.fasta",
		sort_bam = expand("data/raw/sorted_reads/{sample}.sorted.bam", sample = SAMPLES),
		bamqc = expand("reports/bamqc/{sample}_stats/qualimapReport.html", sample = SAMPLES),
		multibamqc = "reports/multisampleBamQcReport.html",
		# CALLING
		hap_caller = expand("data/interim/gvcf_files_bpres/{sample}.raw.snps.indels.g.vcf", sample = SAMPLES),
		split_intervals = expand("data/processed/scattered_intervals/{intervals}-scattered.intervals",
		intervals = INTERVALS),
		joint_geno = expand("data/raw/vcf_bpres/{intervals}.raw.snps.indels.vcf", intervals = INTERVALS),
		# FILTERING
		filtered_snps = expand("data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf", chr = CHR),
		filtered_invariant = expand("data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.invariant.sites.vcf", chr = CHR),
		merge_filtered_snps = "data/processed/filtered_vcf_bpres/allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz",
		merge_filtered_allsites = expand("data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.allsites.vcf.gz", chr = CHR),
		ind_stats = "reports/filtered.qual.dp5_200.maxnocall10.imiss",
		# POPULATION STRUCTURE
		admix_input = "data/processed/biallelic_snps.geno10.maf05.ldpruned.bed",
		admixture = expand("models/admixture/biallelic_snps.geno10.maf05.ldpruned.{k}.Q", k = range(2, 15)),
		pixy_pi = expand("models/pixy/B_oleracea_grouped_{chr}_pi.txt", chr = CHR),
		# DEMOGRAPHY
		vcf2smc = smc_input_files
		# window_pi = expand("models/nucleotide_diversity/{chr}_{population}.windowed.pi", chr = CHR, population = POP)
		# diagnostics = expand("reports/filtering/gvcf_{chr}.table", chr = chr),
		# filter_nocall = expand("data/processed/filtered_snps_bpres/{chr}.filtered.nocall.vcf", chr = chr),
		# diagnostics2 = expand("reports/filtering_bpres/gvcf_{chr}.filtered", chr = chr),
		# merge_vcfs = "data/processed/filtered_snps_bpres/oleracea_filtered.vcf.gz"
		# vcf2smc = expand("models/smc/input/{pop}.{chr}.smc.gz", pop = pops, chr = chr),
		# cv = expand("models/smc/cv/{pop}/fold{fold}/model.final.json", pop = pops, fold = ['0','1']),
		# estimate = expand("models/smc/estimate/{pop}/model.final.json", pop = pops),
		# admix_input = "models/admixture/combined.pruned.bed",
		# admixture = expand("models/admixture/biallelic_snps_geno0.1.pruned.{k}.Q", k = range(0, 15))
		# depth = expand("reports/filtering/dp_{chr}.filtered.dp6_200", chr = chr),
		# filter_depth = expand("data/processed/filtered_snps/{chr}.filtered.dp6_200.nocall.snps.vcf", chr = chr)


################################################################################
## Rule files to include
################################################################################

include: "rules/01-fastqc.smk"
include: "rules/02-mapping.smk"
include: "rules/03-calling.smk"
include: "rules/04-filtering.smk"
include: "rules/05-pop_structure.smk"
include: "rules/06-demography.smk"


### intermediate steps (need to clean up & only include essential outputs in rule_all)

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
