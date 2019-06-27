## SNP calling pipeline for Brassica oleracea
## Follows GATK best practices https://software.broadinstitute.org/gatk/documentation/article.php?id=3893
## Sarah Turner-Hissong
## 4 March 2019

# requires a conda environment to run:
# conda env create --name bo-demography --file environment.yaml

# to run workflow use submit.sh script

rundate = '20190627'
dataset = "cheng2016"

# get sample names (option 1: pull unique identifier from fastq files)
# SAMPLES = glob_wildcards("data/external/fastq_raw/{sample}_pass_1.fastq.gz").sample

# get sample names (option 2: manual entries)
SAMPLES = ['SamC_' + str(x).rjust(3, '0') for x in range(1,120)]
SAMPLES2 = SAMPLES + ['B_rapa']

# dictionary of SRA identifiers from downloaded project info
import csv

with open('data/external/Sra_oleracea.csv', mode = 'r') as infile:
	reader = csv.reader(infile)
	sample_dict = {rows[29]:rows[0] for rows in reader}

# add Brassica rapa reference (Chiifu)
sample_dict.update({'B_rapa':'SRR7881031'})

# list chromosome ids
chr = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']

# add paths for locally installed software
# ngsTools (https://github.com/mfumagalli/ngsTools)
ngsTools = "../software/ngsTools"
# FINESTRUCTURE software (v4)
# wget https://people.maths.bris.ac.uk/~madjl/finestructure/fs_4.0.1.zip
# unzip fs_4.0.1.zip
fs = "../bin/fs_4.0.1"
snpable = "../software/seqbility-20091110"

# create all of the dictionaries
# sample ids for each population
capitata = ','.join(['SamC_' + str(x).rjust(3, '0') for x in range(1,46)])
gongylodes = ','.join(['SamC_' + str(x).rjust(3, '0') for x in range(46,65)])
italica = ','.join(['SamC_' + str(x).rjust(3, '0') for x in range(85,108)])
botrytis = ','.join(['SamC_' + str(x).rjust(3, '0') for x in range(65,85)])
alboglabra = ','.join(['SamC_108', 'SamC_109', 'SamC_110', 'SamC_111'])
gemmifera = ','.join(['SamC_112', 'SamC_113'])
acephala = ','.join(['SamC_114', 'SamC_115'])
sabellica = ','.join(['SamC_116', 'SamC_117'])
wild = ','.join(['SamC_118', 'SamC_119'])

# dictionary of population and corresponding samples
popdict = {
'capitata':'capitata:' + capitata,
'gongylodes':'gongylodes:' + gongylodes,
'italica':'italica:' + italica,
'botrytis':'botrytis:' + botrytis,
'alboglabra' : 'alboglabra:' + alboglabra,
'gemmifera' : 'gemmifera:' + gemmifera,
'acephala' : 'acephala:' + acephala,
'sabellica' : 'sabellica:' + sabellica,
'wild' : 'wild:' + wild,
}

# dictionary of population pairs to estimate divergence
pop_pair_dict = {'botrytis_italica':['botrytis', 'italica'],
'italica_botrytis':['italica', 'botrytis'],
'italica_wild':['italica', 'wild'],
'wild_italica':['wild', 'italica']}

def pop_choose(WC):
	list = popdict[WC.pop]
	return list

# exclude thinning parameters for now
# --thinning controls frequency that CSFS is emitted
# decreasing can potentially increase accuracy of results in recent past
# Default is 1000*log(n), change it to 100*log(n) to provide better resolution in more recent past
# n = haploid sample size in undistinguished portion
# use 100*log((N*2)-1)
# thinning_dict = {
# "capitata": ['195'],
# "gongylodes": ['157'],
# "italica": ['165'],
# "botrytis": ['159'],
# "alboglabra": ['117'],
# "gemmifera": ['60'],
# "acephala": ['60'],
# "sabellica": ['60'],
# "wild": ['60'],
# "kale": ['120']}

# def thinning_choose(WC):
# 	thinning_k = thinning_dict[WC.pop]
# 	return thinning_k

# for split time:
# based on current wildcard, find pop pair in pop_pair_dict
# for each population in current pop pair, find pop in popdict
# combine strings to create a super long string
# return the super long string

def pair_string_choose(WC):
	pops = pop_pair_dict[WC.pop_pair]
	pop1 = popdict[pops[0]]
	pop2 = popdict[pops[1]]
	pop_pair_string = pop1 + " " + pop2
	return pop_pair_string

# create string: location of model.final.json for each pop

def model_chooser(WC):
	pops = pop_pair_dict[WC.pop_pair]
	pop1 = pops[0]
	pop2 = pops[1]
	model_string = "models/estimate/" + WC.rundate + "/" + pop1 + "/model.final.json " \
	+ "models/estimate/" + WC.rundate + "/" + pop2 + "/model.final.json"
	return model_string

# a pseudo-rule that collects the target files (expected outputs)
rule all:
	input:
		# get_ref
		expand("data/external/ref/Boleracea_chromosomes.fasta"),
		# download SRA data, trim ends, and map to reference
		# expand("data/interim/mapped_reads/{sample}.bam", sample = SAMPLES2),
		# sort_bam - output is sorted BAM files for each sample
		expand("data/raw/sorted_reads/{sample}.sorted.bam", sample = SAMPLES2),
		# mark_dups - output is BAM files with duplicates marked
		# expand("data/interim/mark_dups/{sample}.dedup.bam", sample = SAMPLES),
		# add_rg - add read group information to BAM files
		# expand("data/interim/add_rg/{sample}.rg.dedup.bam", sample = SAMPLES),
		# coverage_depth
		expand("reports/bamqc/{sample}_stats/qualimapReport.html", sample = SAMPLES2),
		"reports/multisampleBamQcReport.html",
		# hap_caller - output is GVCF file for each sample
		expand("data/interim/gvcf_files/{sample}.raw.snps.indels.g.vcf", sample = SAMPLES),
		# combine_gvcfs - database for combining multiple gvcf files
		# expand("data/interim/combined_database/{chr}/vcfheader.vcf", chr = chr),
		# directory(expand("data/interim/combined_database/{chr}", chr = chr))
		# joint_geno - outputs joint SNP calls for gvcf files
		expand("data/raw/{chr}.raw.snps.indels.vcf", chr = chr),
		# get_snps
		# expand("data/processed/{chr}.filtered.snps.vcf", chr = chr),
		# bgzip_vcf
		expand("data/processed/{chr}.filtered.snps.vcf.gz", chr = chr),
		# merge vcfs
		"data/processed/merged.vcf.gz",
		# make_kmers
		# expand("data/processed/masks/Boleracea_{chr}.mask.bed.gz", chr = chr),
		# smc
		expand("models/smc/input/{rundate}.{dataset}.{pop}.{chr}.smc.gz",
		rundate = rundate, dataset = 'cheng2016', pop = ['capitata'], chr = chr),
		expand("models/smc/cv/{rundate}.{dataset}.{pop}/model.final.json",
		rundate = rundate, dataset = 'cheng2016', pop = ['capitata'], chr = chr),
		expand("reports/figures/smc/{rundate}.{dataset}.{pop}.png",
		rundate = rundate,
		dataset = 'cheng2016',
		pop = ['capitata'])
		# admix_input
		# "models/admixture/combined.pruned.bed",
		# admixture
		# expand("models/admixture/combined.pruned.{k}.Q", k = [1,2,3,4,5,7,8])
		# # angsd_depth
		# expand("reports/ALL.{chr}.qc.depthGlobal", chr = chr)


include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/filtering.smk"
include: "rules/mask.smk"
include: "rules/demography.smk"
# include: "rules/admixture.smk"
# include: "rules/angsd.smk"
