# run FASTQC for quality control checks on raw sequence data
# exports an html report and gzipped folder with graphs/tables
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

rule fastqc:
	input:
		"data/raw/fastq/{sample}.fastq.gz"
	output:
		"qc/fastqc/{sample}_fastqc.html",
		"qc/fastqc/{sample}_fastqc.zip"
	params:
		"qc/fastqc/"
	threads: 4
	run:
		shell("fastqc \
		-t {threads} \
		-o {params} \
		--noextract \
		{input}")

# aggregate FASTQC reports with multiqc
# generates a html report with quality statistics for all samples (it's awesome)
# chose to generate a separate report for each lane of data,
# but could also generate a single report for all samples
# https://multiqc.info/

rule multiqc:
	input:
		expand("qc/fastqc/{sample}_fastqc.zip", sample = FASTQ)
	output:
		"qc/STJRI01_multiqc.html",
		"qc/STJRI02_multiqc.html",
		"qc/STJRI03_multiqc.html"
	params:
		indir = "qc/fastqc",
		outdir = "qc/"
	run:
		shell("multiqc {params.indir}/*L004* -o {params.outdir} -n STJRI01_multiqc.html")
		shell("multiqc {params.indir}/*L001* -o {params.outdir} -n STJRI02_multiqc.html")
		shell("multiqc {params.indir}/*L002* -o {params.outdir} -n STJRI03_multiqc.html")
