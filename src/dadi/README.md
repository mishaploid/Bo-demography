# Demographic Analyses with dadi

This folder contains the code needed to run demographic inference on the different subpopulations
of *B. oleracea* and its wild relatives.

The main file containing information on the demographic models to run is the `model_config.json` file.
It contains sampling information for each of the *B. oleracea* and wild relative populations, as well as
configuration information for each of the demographic models that will be run.

## Generating site frequency spectra

The `build_frequency_spectra.py` script will generate a dadi-formatted SFS file in the
`sfs/` subfolder for each model specified in `model_config.json`. This script should only need to be run once,
and it will generate all of the SFS files. The input for the script is a VCF file and a populaiton information file,
which is a two-column file containing the names of the samples in the first column and the name of the population
to which they belong in the second column. The name of the output SFS file is specified in `model_config.json`.

**Usage**

```bash
python3 build_frequency_spectra.py \
    --vcf <path-to-vcf-file> \
    --pop_info <path-to-popinfo-file> \
    --subsample
```

## Demographic models

Demographic models are specified in the `demographic_models.py` file.
Each model is a Python function that follows a similar template, taking
in demographic model parameters, the sample sizes of each population, and
the set of grid points to be used for generating the distribution of allele
frequencies. The output of the model function is a SFS, and the demographic history
used to generate the SFS is specified using functions from dadi.

Aside from the demographic models themselves, this script also stores each function in a
dictionary that is used in the `run_inference.py` script to lookup the correct function
when passing in the model name via the command line (see subsection below for more details).

## Inference

The `run_inference.py` script will perform parameter inference using dadi's likelihood framework
for a given demographic model. The model can be specified at the command line using the model's name
as specified in the `model_config.json` file. The model's name is also what will be used to look
up the corresponding demographic model function in the `demographic_models.py` file using the `models`
dictionary.

All of the necessary information for running parameter inference, including the number of likelihood searches to conduct
(default is 100), as well as the initial parameter values and their upper and lower bounds are all pulled from
`model_config.json`.

**Usage**

```bash
python3 run_inference.py --model "cap_gem_vir"
```