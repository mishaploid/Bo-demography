# code to generate indices for reference genome
# source: ftp://ftp.ensemblgenomes.org/pub/plants/release-41/fasta/brassica_oleracea/dna/

# bwa index
# -a bwtsw specifies algorithm for larger genomes
bwa index -a bwtsw Brassica_oleracea.v2.1.dna.toplevel.fa

# fasta indexing
samtools faidx Brassica_oleracea.v2.1.dna.toplevel.fa

# sequence library dictionary
java -jar ~/data/software/picard.jar CreateSequenceDictionary \
    REFERENCE=Brassica_oleracea.v2.1.dna.toplevel.fa \
    OUTPUT=Brassica_oleracea.v2.1.dna.toplevel.dict
