Bo_demography
==============================

Demographic analysis of Brassica oleracea

Genome Alignment 
------------
Follows GATK4 best practices (https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145) 

Implemented as a Snakemake workflow (https://snakemake.readthedocs.io/en/stable/index.html)  

#### STATUS
Works up to HaplotypeCaller step  
Need to manually download FASTQ files, then unhash trimmomatic_pe & bwa_map rules in rules/mapping.smk

#### TODO
* Update mapping.smk (fastq-dump, trimmomatic, and bwa mem)
    + use temp() to set inputs/outputs as temporary so they are removed after running (snakemake currently views as missing inputs)
    + add shell command for fastq-dump 
* Add a config file to specify samples/reference/paths/etc. (will make workflow more general)
* Troubleshoot/optimize calling.smk (works up to HaplotypeCaller step)

### Workflow 

**Snakefile** provides general organization (names of samples, calls other scripts/rules)

**/rules** contains rules (*.smk) that tell snakemake what jobs to run; organized broadly into categories

1. **_mapping.smk_** - map reads to reference genome 
    + trim adapters (Trimmomatic PE)
    + map to reference genome (bwa mem; B. oleracea TO1000 reference)
    + sort BAM files (gatk/SortSam) 
    + mark duplicates with Picard (gatk/MarkDuplicates) 
    + add read groups with Picard (gatk/AddOrReplaceReadGroups)
    
2. **_calling.smk_** - call variants (in progress)
    + identify raw SNPs/indels (gatk/HaplotypeCaller)
    + combine GVCFs (gatk/GenomicsDBImport)
    + joing genotyping (gatk/GenotypeGVCFs)
    
3. **_filtering.smk_** - (in progress) - 
    + select SNPs (gatk/SelectVariants)
    + apply SNP filters (gatk/VariantFiltration)
    
**To run the workflow:**

1. Download FASTQ files to `data/external/fastq_raw` using `fastq-dump`

2. Create a conda environment (only need to do this once)  
`conda env create bo-demography --file environment.yaml`
    + to access for future use, start a screen session with `screen` and use:  
    `source activate bo-demography`

3. If running on a cluster, update the cluster config file (submit.json) 

4. Dry run to test that everything is in place... 
`snakemake -n`

5. Submit the workflow with `./submit.sh` 



Project Organization
------------

    ├── LICENSE
    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── models             <- Trained and serialized models, model predictions, or model summaries
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    ├── src                <- Source code for use in this project.
    │   ├── data           <- Scripts to download or generate data
    │   ├── features       <- Scripts to turn raw data into features for modeling
    │   ├── models         <- Scripts to train models and then use trained models to make
    │   │                     predictions
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations

--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
