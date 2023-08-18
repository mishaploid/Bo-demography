Bo_demography
==============================

Demographic analysis of Brassica oleracea

Alignment of sequence data
------------
Follows GATK4 best practices (https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145)

Implemented as a Snakemake workflow (https://snakemake.readthedocs.io/en/stable/index.html)  

#### TODO
* Clean up filtering steps 
* Add adjustable parameters & metrics for filtering criteria (depth, missingness)  
* More efficient use of scratch space  
* Include sample ids & read groups in bwa-mem step (thanks, Michelle!): bwa mem -R $(echo "@RG\tID:${LANE}_${SAMPLE}\tSM:$SAMPLE\tLB:${SAMPLE}.1\tPL:ILLUMINA") -a $B73FA $R1 $R2  
* Add a config file to specify samples/reference/paths/etc. (will make workflow more general)
  - include variable for scratch directory
* add rule for generating mappability mask using SNPable
* Find a better solution to create input for ADMIXTURE (right now requires manual edit of chr names and snp ids)

### Workflow

**Snakefile** provides general organization (names of samples, calls other scripts/rules)

**/rules** contains rules (*.smk) that tell snakemake what jobs to run; organized broadly into categories

1. **_mapping.smk_** - map reads to reference genome
    + download fastq files (`fasterq-dump`)
    + trim adapters (`Trimmomatic PE`)
    + map to reference genome (`bwa mem`; B. oleracea HDEM reference)
    + sort BAM files (`gatk/SortSam`)
    + mark duplicates with Picard (`gatk/MarkDuplicates`)
    + add read groups with Picard (`gatk/AddOrReplaceReadGroups`)

2. **_calling.smk_** - call variants
    + identify raw SNPs/indels (`gatk/HaplotypeCaller`)
    + combine GVCFs (`gatk/GenomicsDBImport`)
    + joing genotyping (`gatk/GenotypeGVCFs`)

3. **_filtering.smk_** - filter SNPs
    + select SNPs (`gatk/SelectVariants`)
    + apply hard filters (`gatk/VariantFiltration`)
    + export diagnostics (qual, quality depth, depth, marker quality, etc.) (`gatk/VariantsToTable`)
    + remove sites with missing data and select biallelic SNPs (`gat/SelectVariants`)
    + filter on individual sample depth (`gatk/VariantFiltration` + `gatk/SelectVariants`) 
    + export filtered depth for each accession (`gatk/VariantsToTable`)
    + combine vcfs into a single file (`bcftools concat`)
    
4. **_smc.smk_** - estimate population size history with `SMC++`
    + convert vcf to smc format - requires mappability mask (`smc++ vcf2smc`)
    + fit population size history with cross-validation (`smc++ cv`)
    + fit population size history without cross-validation (`smc++ estimate`)
    
5. **_admixture.smk_** - maximum likelihood estimation of individual ancestries
    + create input file (some bash commands + LD pruning in `plink`)
    + run admixture (`admixture`)

**To run the workflow:**

1. Create a conda environment (only need to do this once)  
`conda env create bo-demography --file environment.yaml`
    + to access for future use, start a screen session with `screen` and use:  
    `source activate bo-demography`

2. If running on a cluster, update the cluster config file (submit.json)

3. Dry run to test that everything is in place...
`snakemake -n`

4. Submit the workflow with `./submit.sh`



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
    ├── src                <- Source code for use in this project.
    │   ├── data           <- Scripts to download or generate data
    │   ├── features       <- Scripts to turn raw data into features for modeling
    │   ├── models         <- Scripts to train models and then use trained models to make
    │   │                     predictions
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations

--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
