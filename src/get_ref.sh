#### Get B. oleracea reference v2.1
#### Accession: TO1000
#### See Parkin et al. 2014 Genome Biology

# download unmasked toplevel genomic sequence data
wget -N ftp://ftp.ensemblgenomes.org/pub/plants/release-40/fasta/brassica_oleracea/dna/Brassica_oleracea.v2.1.dna.toplevel.fa.gz

gunzip Brassica_oleracea.v2.1.dna.toplevel.fa.gz
bgzip Brassica_oleracea.v2.1.dna.toplevel.fa

# generate the BWA index (collection of files used by BWA for alignment) 
bwa index -a bwtsw Brassica_oleracea.v2.1.dna.toplevel.fa.gz

# generate fasta file index (output: reference.fa.fai)
# one record per line for each contig (includes contig name, size, location, basesPerLine, and bytesPerLine)
samtools faidx Brassica_oleracea.v2.1.dna.toplevel.fa.gz

# generate sequence dictionary
java -jar ../../../software/picard.jar CreateSequenceDictionary \
    REFERENCE=Brassica_oleracea.v2.1.dna.toplevel.fa.gz \
    OUTPUT=Brassica_oleracea.v2.1.dna.toplevel.dict

