# Selective sweep detection with RAiSD
# https://www.nature.com/articles/s42003-018-0085-8
# https://github.com/alachins/raisd

################################################################################
# STEP 1: generate a bed file of sites that pass filtering criteria

# individual population pruned sample lists generated in R:
#   library(tidyverse)
#   samps <- read_delim("../../pruned_sample_ids.txt", delim = "\t", col_names = c("sample_id", "population"))
#   samps %>% group_by(population) %>% group_walk(~write.table(.x, paste0(.y$population, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE))

# making bed file - need bed file of chromosome lengths:
# cut -f1,1,3 ../data/external/ref/Boleracea_intervals.bed > ../data/external/ref/Boleracea_chr_lengths.bed
# get sites included in vcf file
# get excluded sites that didn't pass filtering criteria/had no coverage

rule create_bed:
    input:
        allsites_vcf = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.allsites.vcf.gz",
        chr_lengths = "data/external/ref/Boleracea_chr_lengths.bed"
    output:
        "models/RAiSD/{chr}_excluded_regions.bed"
    params:
        chr = "{chr}"
    shell:
        """
        zcat {input.allsites_vcf} | sed -e 's/chr//' | awk '{{OFS="\\t"; if (!/^#/){{print $1,$2,$2}}}}' | bedtools complement -i stdin -g {input.chr_lengths} | grep {params.chr} > {output}
        """

# note: to use gzipped vcf as input, RAiSD must be installed with the install-RAiSD-ZLIB.sh script
# -R = include additional information in the report file(s). i.e., window start and end, mu-stat factors
# -s = generate a separate report file per set
# -m = threshold to exclude snps with minor allele frequency < threshold
# -f = force overwrite of previous results with same run ID 
# -X = path to tab-delimited file that contains regions per chromosome to be excluded
# -n = name of the runid
# -I = input VCF file
# -S = sample file, tab delimited 
# -w = window size in SNPs (integer). Default is 50. Only even numbers can be used. 
# -COT = provides a cut-off threshold for identifying top outliers per report. 
# -M = missing data handling strategy; 3 = create a mask for valid alleles and ignore allele pairs with N.
# -y = ploidy level (needed when using -M flag)
rule raisd:
    input:
        allsites_vcf = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.allsites.vcf.gz",
        samples = "models/RAiSD/pruned_sample_lists/{population}.txt",
        excluded_sites = "models/RAiSD/{chr}_excluded_regions.bed" 
    output:
        "models/RAiSD/RAiSD_Report.{population}_{chr}_{window_size}bp.{chr}"
    params:
        runid = "{population}_{chr}_{window_size}bp",
        chr = "{chr}",
        window_size = "{window_size}", # in basepairs
        outdir = "models/RAiSD/"
    resources:
        mem_mb = 50000
    run:
        shell("RAiSD -R -s -f -m 0.05 -X {input.excluded_sites} \
        -n {params.runid} \
        -I {input.allsites_vcf} \
        -S {input.samples} \
        -w {params.window_size} \
        -COT 0.05 \
        -M 3 \
        -y 2")
        shell("mv RAiSD*{params.runid}* {params.outdir}")

# are results comparable across populations if only including doubletons in some? 
# 	-c	Provides the slack for the SFS edges to be used for the calculation of mu_SFS. The default value is 1 (singletons and S-1 snp class, where S is the sample size).

rule raisd_doubletons:
    input:
        allsites_vcf = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.allsites.vcf.gz",
        samples = "models/RAiSD/pruned_sample_lists/{population}.txt",
        excluded_sites = "models/RAiSD/{chr}_excluded_regions.bed" 
    output:
        "models/RAiSD_w_doubletons/RAiSD_Report.{population}_{chr}_{window_size}bp_c2.{chr}"
    params:
        runid = "{population}_{chr}_{window_size}bp_c2",
        chr = "{chr}",
        window_size = "{window_size}", # in basepairs
        outdir = "models/RAiSD_w_doubletons/"
    resources:
        mem_mb = 50000
    run:
        shell("RAiSD -R -s -f -m 0.05 -X {input.excluded_sites} \
        -n {params.runid} \
        -I {input.allsites_vcf} \
        -S {input.samples} \
        -w {params.window_size} \
        -c 2 \
        -COT 0.05 \
        -M 3 \
        -y 2")
        shell("mv RAiSD*{params.runid}* {params.outdir}")
    