# automate download of SRA archive
# www.sthda.com/english/wiki/print.php?id=51

# use command line arguments
args <- commandArgs(trailingOnly = TRUE)

# read in SRA info
sri <- read.csv("data/external/Sra_oleracea.csv", stringsAsFactors = FALSE)
files <- basename(sri$download_path)

# extract accession names
# accession <- gsub("Brassica oleracea ", "", sri[,30])
accession <- args[1]
SRR <- sri$Run[sri$SampleName == paste0("Brassica oleracea ", accession)]

print(accession)
print(SRR)

# send command to terminal
cmd = paste("fasterq-dump", SRR, "-O data/external/fastq_raw -o", accession, "-t /scratch/sdturner -p")
# see fastq-dump -h for details

cat(cmd, "\n") # print current command

system(cmd) # invoke command
