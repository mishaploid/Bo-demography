
# run ADMIXTURE analysis

# format input data
# combine chromosome vcfs to single vcf
# convert to bed/bim/fam format
# LD pruning

# # phasing - beagle 4.1
# # consider adjusting window size and ibd paramaters
# rule phase_vcf:
#     input:
#         vcf = "data/processed/{chr}.filtered.snps.vcf.gz"
#         map =
#     output:
#         "models/beagle/{chr}.phased.filtered.snps.vcf.gz"
#     params:
#         out_prefix = "models/beagle/{chr}.phased.filtered.snps"
#         chrom = "{chr}"
#         ibd = "true"
#     threads: 32
#     run:
#         shell("beagle \
#         nthreads={threads}\
#         gt={input.vcf} \
#         chrom={params.chrom} \
#         out={params.out_prefix} \
#         map={input.map} \
#         ibd={params.ibd}")

# need to use plink2 (see http://apol1.blogspot.com/2014/11/best-practice-for-converting-vcf-files.html)
# download:
# cd ~/bin
# wget http://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20190402.zip
# unzip plink2_linux_x86_64_20190402.zip
# add `export PATH="/home/sdturner/bin:$PATH"` to ~/.bashrc profile

# SCRATCH
# ~/software/plink2 --vcf ../../data/processed/filtered_snps_bpres/all_samps.filtered.snps.invariant.vcf.gz --allow-extra-chr --max-alleles 2 --mind 0.2 --geno 0.2 --make-bed --out filtered.mind.2.geno.2

rule admix_input:
    input:
    	ref = "data/external/ref/Boleracea_chromosomes.fasta",
        vcf = expand("data/processed/filtered_snps/{chr}.filtered.dp6_200.nocall.snps.vcf.gz", chr = chr)
    output:
    	"models/admixture/biallelic_snps_geno0.1.pruned.bed"
    params:
        vcf = "data/processed/filtered_snps/oleracea_filtered.vcf.gz",
        stem = "models/admixture/biallelic_snps_geno0.1",
        pruned = "models/admixture/biallelic_snps_geno0.1.pruned"
    run:
		shell("bcftools concat {input.vcf} -Oz -o {output}")
		shell("tabix -p vcf {output}")
        shell("~/software/plink2 --vcf {params.vcf} \
        --allow-extra-chr \
        --max-alleles 2 \
        --vcf-filter \
        --geno 0.1 \
        --make-bed \
        --out {params.stem}")
        shell("~/software/plink2 --bfile {params.stem} \
        --indep-pairwise 50 10 0.1 \
        --out {params.stem}")
        shell("plink2 --bfile {params.stem} \
        --extract {params.stem}.prune.in \
        --out {params.pruned} \
        --make-bed")
        shell("plink --bfile {params.stem} \
        --extract {params.stem}.prune.in \
        --out {params.pruned} \
        --recode")

# could use distruct for plotting...

rule admixture:
    input:
        bed = "models/admixture/biallelic_snps_geno0.1.pruned.bed"
    output:
        "models/admixture/biallelic_snps_geno0.1.pruned.{k}.Q",
        "models/admixture/biallelic_snps_geno0.1.pruned.{k}.P"
    params:
        k = "{k}"
    threads: 32
    run:
        shell("admixture -B --cv -j{threads} {input.bed} {params.k}")
        shell("mv combined.pruned.{params.k}.* models/admixture")
