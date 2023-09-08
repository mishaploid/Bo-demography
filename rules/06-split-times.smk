# # #Generate vcf2smc files containing the joint frequency spectrum for both populations

rule joint_vcf2smc12:
    input:
        vcf = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz",
        index = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz.tbi",
        mask = "data/processed/mappability_masks/scratch/{chr}.mask.bed.gz"
    output:
        out12 = "models/smc_split_input/{pop_pair}_12.{distinguished_ind1}.{chr}.smc.gz",
    params:
    	chrom = "{chr}",
        distind1 = "{distinguished_ind1}",
    	pop_pair_string12 = pair_string_choose12
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
        """
        smc++ vcf2smc \
        --mask {input.mask} \
        -d {params.distind1} {params.distind1} \
        {input.vcf} {output.out12} {params.chrom} {params.pop_pair_string12}
        """

rule joint_vcf2smc21:
    input:
        vcf = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz",
        index = "data/processed/filtered_vcf_bpres/{chr}_allsamps.filtered.qual.dp5_200.maxnocall10.biallelic.snps.vcf.gz.tbi",
        mask = "data/processed/mappability_masks/scratch/{chr}.mask.bed.gz"
    output:
        out21 = "models/smc_split_input/{pop_pair}_21.{distinguished_ind2}.{chr}.smc.gz"
    params:
    	chrom = "{chr}",
        distind2 = "{distinguished_ind2}",
        pop_pair_string21 = pair_string_choose21
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
        """
        smc++ vcf2smc \
        --mask {input.mask} \
        -d {params.distind2} {params.distind2} \
        {input.vcf} {output.out21} {params.chrom} {params.pop_pair_string21}
        """

rule smc_split:
    input:
        smc_split_input
    params:
        model_out_dir = "models/smc_split/{pop_pair}/"
    output:
        model_out = "models/smc_split/{pop_pair}/model.final.json"
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
        "smc++ split \
        --cores 8 \
        -o {params.model_out_dir} \
        {input}"


# bootstrapping
rule joint_bootstrap_vcf2smc12:
    input:
        expand("models/smc_split_input/{{pop_pair}}_12.{{distinguished_ind1}}.{chr}.smc.gz", chr = CHR)
    output:
        expand('models/smc_split_bootstrap_input/{{pop_pair}}_12.{{distinguished_ind1}}_rep_{n_bootstrap}/bootstrap_chr{boot_chr}.gz', n_bootstrap = range(1,11), boot_chr = range(1,10)),
    params:
    	pop_pair = "{pop_pair}",
        distind1 = "{distinguished_ind1}",
        nr_bootstraps = 10,
        chunk_size = 5000000,
        nr_chr = 9,
        input_dir = "models/smc_split_input/{pop_pair}_12.{distinguished_ind1}*"
    shell:
        "python3 scripts/smc_bootstrap.py \
        --nr_bootstraps {params.nr_bootstraps} \
        --chunk_size {params.chunk_size} \
        --chunks_per_chromosome 10 \
        --nr_chromosomes {params.nr_chr} \
        models/smc_split_bootstrap_input/{params.pop_pair}_12.{params.distind1}_rep \
        {params.input_dir}"

rule joint_bootstrap_vcf2smc21:
    input:
        expand("models/smc_split_input/{{pop_pair}}_21.{{distinguished_ind1}}.{chr}.smc.gz", chr = CHR)
    output:
        expand('models/smc_split_bootstrap_input/{{pop_pair}}_21.{{distinguished_ind1}}_rep_{n_bootstrap}/bootstrap_chr{boot_chr}.gz', n_bootstrap = range(1,11), boot_chr = range(1,10)),
    params:
    	pop_pair = "{pop_pair}",
        distind1 = "{distinguished_ind1}",
        nr_bootstraps = 10,
        chunk_size = 5000000,
        nr_chr = 9,
        input_dir = "models/smc_split_input/{pop_pair}_21.{distinguished_ind1}*"
    shell:
        "python3 scripts/smc_bootstrap.py \
        --nr_bootstraps {params.nr_bootstraps} \
        --chunk_size {params.chunk_size} \
        --chunks_per_chromosome 10 \
        --nr_chromosomes {params.nr_chr} \
        models/smc_split_bootstrap_input/{params.pop_pair}_21.{params.distind1}_rep \
        {params.input_dir}"

rule smc_split_bootstrap:
    input:
        smc_split_bootstrap_input
    params:
        model_out_dir = "models/smc_split_bootstrap/{pop_pair}_{n_bootstrap}"
    output:
        model_out = "models/smc_split_bootstrap/{pop_pair}_{n_bootstrap}/model.final.json"
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
        "smc++ split \
        --cores 2 \
        -o {params.model_out_dir} \
        {input}"

rule plot_split:
    input:
        smc_split = expand("models/smc_split/{pop_pair}/model.final.json", pop_pair = pop_pair_dict.keys()),
        smc_split_bootstrap = expand("models/smc_split_bootstrap/{pop_pair}_{n_bootstrap}/model.final.json", pop_pair = pop_pair_dict.keys(), n_bootstrap = range(1,11))
    output:
        split = "reports/carrot_all_pops_smc_split.png",
        bootstrap = "reports/carrot_all_pops_smc_split_bootstrap.png"
    params:
        gen = config['gen']
    singularity:
        "docker://terhorst/smcpp:latest"
    shell:
        """
        smc++ plot \
        --csv \
        -g {params.gen} \
        {output.split} \
        {input.smc_split}
        smc++ plot \
        --csv \
        -g {params.gen} \
        {output.bootstrap} \
        {input.smc_split_bootstrap}
        """
