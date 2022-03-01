################################################################################
# run FASTQC for quality control checks on raw sequence data
# 	exports an html report and gzipped folder with graphs/tables
# 	https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
################################################################################

rule fastqc:
	input:
		"data/raw/fastq/{sample}_{readgroup}_001.fastq.gz"
	output:
		touch("qc/fastqc/{sample}_{readgroup}_fastqc.html"),
		touch("qc/fastqc/{sample}_{readgroup}_fastqc.zip")
	params:
		"qc/fastqc/"
	threads: 4
	run:
		shell("fastqc \
		-t {threads} \
		-o {params} \
		--noextract \
		{input}")

################################################################################
# aggregate FASTQC reports with multiqc
# 	generates a html report with quality statistics for all samples (it's awesome)
# 	chose here to generate a separate report for each lane of data,
# 	but could also generate a single report for all samples
# 	https://multiqc.info/
################################################################################

rule multiqc:
	input:
		expand("qc/fastqc/{sample}_fastqc.zip", sample = SAMPLES)
	output:
		touch("qc/STJRI01_multiqc.html"),
		touch("qc/STJRI02_multiqc.html"),
		touch("qc/STJRI03_multiqc.html")
	params:
		indir = "qc/fastqc",
		outdir = "qc/"
	run:
		shell("multiqc {params.indir}/*L004* -o {params.outdir} -n STJRI01_multiqc.html")
		shell("multiqc {params.indir}/*L001* -o {params.outdir} -n STJRI02_multiqc.html")
		shell("multiqc {params.indir}/*L002* -o {params.outdir} -n STJRI03_multiqc.html")
