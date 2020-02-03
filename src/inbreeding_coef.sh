#!/bin/bash

set -e
set +o pipefail

## adapted from ANGSD-wrapper
## https://github.com/ANGSD-wrapper/angsd-wrapper/blob/master/Wrappers/Inbreeding_Coefficients.sh

# set directories for ANGSD and NGStools
ANGSD_DIR="../software/angsd"
NGSF_DIR="../software/ngstools/ngsF"

# use inputs from snakemake
INPUT=$(realpath $1)
REGIONS=$(realpath $2)
OUT=$(realpath $3)
REF=$(realpath $4)
ANC=$(realpath $5)

N_IND=$(wc -l < ${INPUT})

N_CORES=32

# get genotype likelihoods from ANGSD
${ANGSD_DIR}/angsd \
    -bam ${INPUT} \
    -rf ${REGIONS} \
    -doGLF 3 \
    -GL 1 \
    -out ${OUT}\
    -ref ${REF} \
    -anc ${ANC} \
    -doMaf 1 \
    -SNP_pval 1e-6 \
    -doMajorMinor 1 \
    -doCounts 1 \
    -minMapQ 30 \
    -minQ 20 \
    -setMinDepth 3 \
    -setMaxDepth 2156 \
    -nThreads ${N_CORES}

# tally number of sites for calculating inbreeding coefficients
N_SITES=`expr $(zcat "${OUT}".mafs.gz | wc -l) - 1`
SEED=$RANDOM # set a random seed

zcat ${OUT}.glf.gz | ${NGSF_DIR}/ngsF \
    -glf - \
    -out ${OUT}.approx_indF \
    -n_ind ${N_IND} \
    -n_sites ${N_SITES} \
    -min_epsilon 1e-9 \
    -approx_EM \
    -seed ${SEED} \
    -init_values r \
    -n_threads ${N_CORES}

zcat ${OUT}.glf.gz | ${NGSF_DIR}/ngsF \
    -glf - \
    -out ${OUT}.indF \
    -n_ind ${N_IND} \
    -n_sites ${N_SITES} \
    -min_epsilon 1e-9 \
    -init_values ${OUT}.approx_indF.pars \
    -n_threads ${N_CORES}
