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

rule raisd:
    input:
        allsites_vcf = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.allsites.vcf.gz",
        samples = "models/RAiSD/pruned_sample_lists/{population}.txt",
        excluded_sites = "models/RAiSD/{chr}_excluded_regions.bed" 
    output:
        "models/RAiSD/RAiSD_Report.{population}_{chr}_w{window_size}kb.{chr}"
    params:
        population = "{population}",
        chr = "{chr}",
        window_size = "{window_size}", # in kb 
        out = "*.{params.runid}_w{window_size}kb.{chr}*",
        outdir = "models/RAiSD/"
    run:
        shell("RAiSD -R -s -m 0.05 -P -f -X {input.excluded_sites} \
        -n {params.runid} \
        -I {input.allsites_vcf} \
        -S {input.samples} \
        -w {params.window_size}")
        shell("mv {params.out}* {params.outdir}")