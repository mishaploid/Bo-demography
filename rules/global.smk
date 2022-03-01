import csv
import os
import pandas as pd

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
pop_ids = pd.read_csv(config['pruned_pop_ids'], sep="\t", names=['Sample_ID', 'Population'])

# filter for pops with more than 4 samples
pop_ids = pop_ids.groupby(['Population']).filter(lambda x: len(x) > 3)

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

# for bootstrapped inputs
smc_bootstrap_input = [expand('models/smc/bootstrap_input/{population}.{distinguished_ind}_rep_{n_bootstrap}/bootstrap_chr{boot_chr}.gz',
                    population = key, distinguished_ind = value, n_bootstrap = range(1,11), boot_chr = range(1,10))
                    for key, value in distind_dict.items()]

# create a list of filenames to feed into smc++ cv command
# want to include all .smc.gz files for each population (includes all chromosomes and distinguished individuals)
def smc_cv_input(wildcards):
    files = expand('models/smc/input/{population}.{distinguished_ind}.{chr}.smc.gz', chr=CHR, population=wildcards.population, distinguished_ind=distind_dict[wildcards.population])
    return files

def smc_cv_boot_input(wildcards):
    files = expand('models/smc/bootstrap_input/{population}.{distinguished_ind}_rep_{n_bootstrap}/bootstrap_chr{boot_chr}.gz', population=wildcards.population, distinguished_ind=distind_dict[wildcards.population], n_bootstrap = range(1,11), boot_chr = range(1,10))
    return files

# population specific sample lists for selective sweeps
# after silhouette pruning
