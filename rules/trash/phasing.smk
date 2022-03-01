### phasing in BEAGLE 4.1
# wanted to use v5 but ran into memory errors
# Note: using default parameters
# ideally would include a linkage map as more precise estimate of recombination

rule beagle_phasing:
    input:
        "data/processed/filtered_snps/{chr}.filtered.snps.vcf.gz"
    output:
        "data/processed/phased/{chr}.phased.filtered.vcf.gz"
    params:
        outstring = "data/processed/phased/{chr}.phased.filtered"
    run:
        shell("java -Xmx16g -jar {beagle} gt={input} out={params.outstring}")
