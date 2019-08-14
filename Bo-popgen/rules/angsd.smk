# STEP 1: assess depth (overall and per individual)

# Note - originally was using ANGSD for coverage, but don't like output format much
# switched to samtools depth & moved to filtering.smk

# STEP 2: create a fasta file for ancestral reference
# input: short reads from B. rapa reference (Chiifu-401) mapped to B. oleracea reference (HDEM)
    # -i                input file (BAM format)
    # -out              prefix for output
    # -ref              reference sequence (B. oleracea HDEM assembly)
    # -dofasta          create fasta file from BAM (option 2: use most common non-missing base)
    # -doCounts         frequency of different bases (must include for -dofasta 1)
    # -uniqueOnly       remove reads that have multiple best hits (0: no, 1: remove)
    # -remove_bads      remove reads with flag above 255 (not primary, failure, duplicate reads)
    # -only_proper_pairs use only reads where both mates mapped correctly (1: include; 0: use all reads)
    # -baq              perform BAQ (base alignment quality) computation (0: none, 1: normal, 2: extended)
    # -minMapQ          filter for minimum mapQ quality (set to 30)
    # -minQ             filter for minimum base quality score (set to 20)

rule anc_fasta:
	input:
		anc = "../data/raw/sorted_reads/B_rapa.sorted.bam",
		ref = "../data/external/ref/Boleracea_chromosomes.fasta"
	output:
		"data/processed/anc/B_rapa_HDEM.fa.gz"
	params:
		out = "data/processed/anc/B_rapa_HDEM"
	shell:
		"angsd \
		-i {input.anc} \
		-out {params.out} \
		-ref {input.ref} \
		-dofasta 2 \
		-doCounts 1 \
		-uniqueOnly 1 \
		-remove_bads 1 \
		-only_proper_pairs 1 \
        -trim 0 \
        -C 50 \
		-baq 1 \
		-minMapQ 30 \
		-minQ 20"

# STEP 3: pcangsd input
# get list of bam files:
# find ../data/raw/sorted_reads |  grep bam$ > data/all.files

rule beaglegl:
	input:
		ref = "../data/external/ref/Boleracea_chromosomes_only.fa",
		bamlist = "data/all.files"
	output:
	    "data/processed/pcangsd_{chr}.beagle.gz"
	params:
		out = "data/processed/pcangsd_{chr}",
		chr = "{chr}"
	threads: 50
	run:
		shell("angsd -GL 1 \
		-out {params.out} \
		-nThreads {threads} \
		-doGlf 2 \
		-doMajorMinor 1 \
		-doMaf 2 \
		-SNP_pval 1e-6 \
		-uniqueOnly 1 \
		-remove_bads 1 \
		-only_proper_pairs 1 \
		-C 50 \
		-baq 1 \
		-minmapQ 30 \
		-minQ 20 \
		-minInd 140 \
		-bam {input.bamlist} \
		-ref {input.ref} \
		-r {params.chr}")
		# shell("python {pcangsd}/pcangsd.py \
		# -beagle {params.out}.beagle.gz \
		# -inbreed 2 \
		# -threads {threads} \
		# -o {params.out2}")


# create input for principal component analysis
#
# rule geno_likelihood:
#     input:
#         ref = "data/external/ref/Boleracea_chromosomes.fasta",
#         bamlist = "data/processed/all.bamlist"
#     output:
#         geno = "reports/PCA/B_oleracea_PCA_{chr}.geno.gz"
#     benchmark:
#         "benchmarks/pca.input.benchmark"
#     params:
#         chr = "{chr}",
#         outfile = "reports/PCA/B_oleracea_PCA_{chr}"
#     shell:
#         "{ANGSD}/angsd \
#         -P 8 \
#         -bam {input.bamlist} \
#         -r {params.chr} \
#         -out {params.outfile} \
#         -ref {input.ref} \
#         -uniqueOnly 1 \
#         -remove_bads 1 \
#         -only_proper_pairs 1 \
#         -GL 2 \
#         -doGLF 1 \
#         -doMaf 1 \
#         -doMajorMinor 1 \
#         -skipTriallelic 1 \
#         -SNP_pval 1e-3 \
#         -minMapQ 30 \
#         -minQ 20 \
#         -nInd 119 \
#         -minInd 100 \
#         -setMinDepth 3 \
#         -setMaxDepth 2156 \
#         -baq 1 \
#         -doCounts 1 \
#         -doGeno 32 \
#         -doPost 1"
#
# rule ngs_pca:
#     input:
#
#     output:
#         covar = "reports/PCA/ALL.covar"
#     shell:
#         "{NGSTOOLS}/ngsPopGen/ngsCovar \
#         -probfile {input.geno} \
#         -outfile {output.covar} \
#         -nind 119 \
#         -nsites 100000 \
#         -call 0 \
#         -norm 0"
#
# rule inb_coef:
# 	input:
# 		bamlist = "data/processed/all.bamlist",
# 		regions = "data/processed/regions.txt",
# 		ref = "data/external/ref/Brassica_oleracea.v2.1.dna.toplevel.fa",
# 		anc = "data/processed/anc/B.rapa.anc.fa"
# 	output:
# 		"reports/Inbreeding_Coefficients/B_oleracea_all.approx_indF"
# 	params:
# 		out = "reports/Inbreeding_Coefficients/B_oleracea_all"
# 	shell:
# 		"""
# 		./src/inbreeding_coef.sh {input.bamlist} {input.regions} {params.out} {input.ref} {input.anc}
# 		"""
#
# # rule inb_coef:
# #     input:
# #         glf = "reports/Inbreeding_Coefficients/B_oleracea_all.glf.gz"
# #     output:
# #         pars = "reports/Inbreeding_Coefficients/B_oleracea_all.approx_IndF.pars",
# #         indF = "reports/Inbreeding_Coefficients/B_oleracea_all.indF"
# #     params:
# #         nsites = "21893327"
# #     run:
# #         shell("zcat {input.glf} | ../software/angsd-wrapper/dependencies/ngsF/ngsF --n_ind 119 --n_sites {params.nsites} --glf - --min_epsilon 0.001 --out B_oleracea_all.approx_indF --approx_EM --seed 12345 --init_values r")
# #         shell("zcat {input.glf} | ../software/angsd-wrapper/dependencies/ngsF/ngsF --n_ind 119 --n_sites {params.nsites} --glf - --min_epsilon 0.001 --out B_oleracea_all.indF --init_values B_oleracea_all.approx_indF.pars")
#
#
# # previous input:
# 		        # bamlist = lambda wildcards: expand('data/processed/{pop}.bamlist',
# 		        #         pop = wildcards.pop),
# 		        # ref = "data/external/ref/Brassica_oleracea.v2.1.dna.toplevel.fa"
# 				# output:         # "reports/neutrality/{pop}.saf.idx"
#
# # note using ref as the ancestral for folded SFS
#
# 		# shell("{ANGSD}/angsd -bam {input.bamlist} -ref {input.ref} -anc {input.ref} -fold 1 \
# 		# -out reports/neutrality/{params.pop} \
# 		# -rf {params.rf} \
# 		# -uniqueOnly 1 \
# 		# -remove_bads 1 \
# 		# -only_proper_pairs 1 \
# 		# -baq 1 \
# 		# -minMapQ 30 \
# 		# -minQ 20 \
# 		# -minInd 1 \
# 		# -setMinDepth 20 \
# 		# -setMaxDepth 2156 \
# 		# -doCounts 1 \
# 		# -GL 1 \
# 		# -doSaf 1")
# 		# shell("{ANGSD}/misc/realSFS reports/neutrality/{params.pop}.saf.idx > {params.sfs}")
#
# rule neutrality_stats:
# 	input:
# 		bamlist = "data/processed/capitata.bamlist",
# 		ref = "data/external/ref/Brassica_oleracea.v2.1.dna.toplevel.fa"
# 	output:
# 		"reports/neutrality/{pop}.{chr}.thetas.idx"
# 	params:
# 		pop = "{pop}",
# 		chr = "{chr}",
# 		rf = "data/processed/regions.txt",
# 		sfs = "reports/neutrality/{pop}.sfs",
# 		tmp = "reports/neutrality/{pop}.{chr}.thetas.idx",
# 		out = "reports/neutrality/{pop}.{chr}"
# 	run:
# 		shell("{ANGSD}/angsd -bam {input.bamlist} -r {params.chr} \
# 		-out {params.out} \
# 		-uniqueOnly 1 \
# 		-remove_bads 1 \
# 		-only_proper_pairs 1 \
# 		-baq 1 \
# 		-doCounts 1 \
# 		-minMapQ 30 \
# 		-minQ 20 \
# 		-setMinDepth 20 \
# 		-setMaxDepth 2156 \
# 		-doThetas 1 \
# 		-doSaf 1 \
# 		-pest {params.sfs} \
# 		-ref {input.ref} \
# 		-anc {input.ref} \
# 		-fold 1 -GL 1")
#
# 		# shell("{ANGSD}/misc/thetaStat do_stat {params.tmp} \
# 		# -nChr 9 \
# 		# -win 50000 \
# 		# -step 10000 \
# 		# -outnames {params.out}")
#
#
# # NOTE: need to update to run .maf and .glf by chromosome first, then combine & run inbreeding coef
# # rule inb_coef:
# #     input:
# #         "src/models/Inbreeding_Coefficients_Config"
# #     output:
# #         indf = "reports/Inbreeding_Coefficients/B_oleracea_all.indF",
# #         pars = "reports/Inbreeding_Coefficients/B_oleracea_all.approx_indF.pars"
# #     log:
# #         "logs/inb_coef.log"
# #     shell:
# #         "{ANGSD} Inbreeding {input}"
#
#     # input:
#     #     ref = "data/external/ref/Brassica_oleracea.v2.1.dna.toplevel.fa",
#     #     bamlist = "data/processed/all.bamlist"
#     # output:
#     #
#     # log:
#     # params:
#     #     chr = "{chr}"
#     # shell:
#     #     "{ANGSD}/angsd -P 8 -b {input.bamlist} -ref {input.ref} "
#     #     "-r {chr} -doGLF 2 -GL 2 "
#     #     "-out {output} "
#     #     "-ref {input.ref} -anc {input.ref} "
#     #     "-doMaf 1 -SNP_pval -doMajorMinor"
#     #     "-minMapQ 30 -minQ 20 -minInd 113 -setMinDepth 3 -setMaxDepth 2156 "
#     #
#
#     # input:
#     #     ref = "data/external/ref/Brassica_oleracea.v2.1.dna.toplevel.fa",
#     #     bamlist = "data/processed/all.bamlist"
#     # output:
#     #     "reports/all.{chr}.qc.arg",
#     #     "reports/all.{chr}.qc.depthGlobal",
#     #     "reports/all.{chr}.qc.depthSample",
#     #     "reports/all.{chr}.qc.qs"
#     # log:
#     #     "logs/qc_filter.log"
#     # params:
#     #     out = "reports/all.{chr}.qc",
#     #     chr = "{chr}"
#     # shell:
#     #     "{ANGSD}/angsd -P 8 -b {input.bamlist} -ref {input.ref} -out {params.out} "
#     #     "-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 "
#     #     "-minMapQ 30 -minQ 20 -minInd 113 -r {params.chr} "
#     #     "-doQsDist 1 -doDepth 1 -doCounts 1"
#
# # Rscript ../../software/ngstools/Scripts/plotQC.R all.C9.qc
