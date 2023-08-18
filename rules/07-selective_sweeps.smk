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
        allsites_vcf = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.allsites.vcf",
        chr_lengths = "data/external/ref/Boleracea_chr_lengths.bed"
    output:
        "models/RAiSD/{chr}_excluded_regions.bed"
    params:
        chr = "{chr}"
    shell:
        """
        cat {input.allsites_vcf} | sed -e 's/chr//' | awk '{{OFS="\\t"; if (!/^#/){{print $1,$2,$2}}}}' | bedtools complement -i stdin -g {input.chr_lengths} | grep {params.chr} > {output}
        """

rule unzip_vcf:
    input:
        "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.allsites.vcf.gz"
    output:
        "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.allsites.vcf"
    run:
        shell("gunzip -c {input} > {output}")

rule raisd:
    input:
        allsites_vcf = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.allsites.vcf",
        samples = "models/RAiSD/pruned_sample_lists/{population}.txt",
        excluded_sites = "models/RAiSD/{chr}_excluded_regions.bed" # "data/processed/mappability_masks/scratch/{chr}.mask.bed.gz"
    output:
        "models/RAiSD/RAiSD_Report.{population}_{chr}_c2_w500.{chr}"
    params:
        population = "{population}_{chr}_c2",
        out = "*.{population}_{chr}_c2_w500.{chr}*",
        outdir = "models/RAiSD/"
    run:
        shell("RAiSD -R -s -m 0.05 -P -X {input.excluded_sites} -n {params.population} -I {input.allsites_vcf} -S {input.samples} -w 500 -f -c 2")
        shell("mv {params.out}* {params.outdir}")


# correct for invariant sites using mop (thanks Silas for the code!!!)
# "mop -m 0.7 -x 100 -i 5 -Q 30 -q 30 -b {input.bams} -R {params.chrom}:{params.start}-{params.end} > {output}"
#
# rule var_correct:
#     input:
#         raisd="data/raisd/RAiSD_Report.{ref}--{ssp}--{pop}--{c}--{r1}--{r2}.txt",
#         bed = "data/mop/{ref}--{ssp}--{pop}--{c}--{r1}--{r2}.bed"
#     output:
#         "data/raisd/RAiSD_Report.{ref}--{ssp}--{pop}--{c}--{r1}--{r2}.corrected"
#     params:
#         c = "{c}"
#     shell:
#         """
#         awk -v c={params.c} 'BEGIN{{OFS = "\\t"}}; NR > 1 {{print c, $2-1, $3, $0, $3-$2+1}}' {input.raisd} |\
#         bedtools groupby -i - -g 1,2,3,4,5,6 -c 7,8,9,10,11 -o last,last,last,last,distinct |\
#         bedtools intersect -a - -b {input.bed} |\
#         awk 'BEGIN{{OFS = "\\t"}};{{print $0, $3-$2}}' |\
#         bedtools groupby -i stdin -g 1,4,5,6 -c 7,8,9,10,11,12 -o distinct,distinct,distinct,distinct,distinct,sum |\
#         awk 'BEGIN{{OFS="\\t"}}{{varnew = ($10/$9)*$5; print $1, $2, $3, $4, varnew, $6, $7, varnew*$6*$7}}' > {output}
#         """
