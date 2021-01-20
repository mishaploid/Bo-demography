
# run ADMIXTURE analysis

# format input data
# convert to bed/bim/fam format (uses plink2)
#   http://apol1.blogspot.com/2014/11/best-practice-for-converting-vcf-files.html
# LD pruning

# SCRATCH
# ~/software/plink2 --vcf ../../data/processed/filtered_snps_bpres/all_samps.filtered.snps.invariant.vcf.gz --allow-extra-chr --max-alleles 2 --mind 0.2 --geno 0.2 --make-bed --out filtered.mind.2.geno.2

rule admix_input:
    input:
    	ref = config['ref'],
        vcf = "data/processed/filtered_vcf_bpres/allsamps.filtered.qual.dp6_200.maxnocall10.biallelic.snps.vcf.gz"
    output:
    	"data/processed/biallelic_snps.geno10.maf05.ldpruned.bed"
    params:
        plink2 = config['plink2'],
        stem = "data/processed/biallelic_snps.geno10.maf05",
        pruned = "data/processed/biallelic_snps.geno10.maf05.ldpruned"
    run:
        shell("{params.plink2} --vcf {input.vcf} \
        --allow-extra-chr \
        --max-alleles 2 \
        --vcf-filter \
        --geno 0.1 \
        --maf 0.05 \
        --make-bed \
        --out {params.stem}")
        shell("cp {params.stem}.bim {params.stem}_temp.bim")
        shell("""sed \"s/^C//g" {params.stem}_temp.bim > {params.stem}_temp2.bim""")
        shell("./src/replace_rsids.sh {params.stem}_temp2.bim > {params.stem}.bim")
        shell("rm {params.stem}_temp*.bim")
        shell("{params.plink2} --bfile {params.stem} \
        --allow-extra-chr \
        --indep-pairwise 50 10 0.1 \
        --out {params.stem}")
        shell("{params.plink2} --bfile {params.stem} \
        --allow-extra-chr \
        --extract {params.stem}.prune.in \
        --out {params.pruned} \
        --make-bed")

# could use distruct for plotting...
# Note default number of bootstraps = 200

rule admixture:
    input:
        bed = "data/processed/biallelic_snps.geno10.maf05.ldpruned.bed"
    output:
        "models/admixture/biallelic_snps.geno10.maf05.ldpruned.{k}.Q",
        "models/admixture/biallelic_snps.geno10.maf05.ldpruned.{k}.P"
    params:
        stem = "biallelic_snps.geno10.maf05.ldpruned",
        k = "{k}"
    threads: 32
    run:
        shell("admixture -B --cv -j{threads} {input.bed} {params.k}")
        shell("mv {params.stem}.{params.k}.* models/admixture")

# calculate nucleotide diversity for 10 kb windows
# in R, to generate sample lists:
#      samps <- read_delim("sample_ids.txt", delim = "\t")
#      samps %>% group_by(pop) %>% group_walk(~write_delim(.x, paste0(.y$pop, "_samps.txt"), delim = "\t", col_names = FALSE))

rule window_pi:
    input:
        allsites_vcf = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp6_200.maxnocall10.allsites.vcf.gz",
        keep = "models/samp_lists/{pop}_samps.txt",
        remove = "models/samp_lists/sil_drop.txt"
    output:
        window_pi = "models/nucleotide_diversity/{chr}_{pop}.windowed.pi"
    params:
        chr = "{chr}"
    run:
        shell("vcftools --gzvcf {input.allsites_vcf} \
        --chr {params.chr} \
        --keep {input.keep} \
        --remove {input.remove} \
        --out {output.window_pi} \
        --window-pi 10000")
