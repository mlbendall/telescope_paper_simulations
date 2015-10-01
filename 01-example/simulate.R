#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

transcripts_fasta <- args[1]
proportion_table  <- args[2]
reads_per_samp    <- as.integer(args[3])
seed              <- as.integer(args[4])
if(is.na(seed)) seed <- NULL

outdir <- 'simdata'

# Load Biostrings
library(Biostrings)

# Load development version of polyester
sourcepath <- paste0( Sys.getenv("POLYESTER_DIR"), '/R')
sourcefiles <- list.files(path=sourcepath,pattern="*.R")
sapply(paste(sourcepath,sourcefiles,sep='/'),source)

# Load the FASTA file
fasta <- readDNAStringSet(transcripts_fasta)

# Load proportions for simulation
props <- read.table(proportion_table, sep='\t', header=T, row.names="locus")

# Subset the FASTA file
small_fasta <- fasta[rownames(props)]
smallseqfile <- tempfile()
writeXStringSet(small_fasta, smallseqfile)

# Sequencing parameters
fraglen        <- 450
fragsd         <- 25
readlen        <- 90
error_rate     <- 0.0
paired         <- 1

cat('reads per sample ', reads_per_samp, '\n')
actual_cov <- reads_per_samp * (paired+1) * readlen / sum(width(small_fasta))
cat('calculated coverage ', actual_cov, '\n')
cat('seed: ', seed, '\n')

experiment_counts <- round(props * reads_per_samp)

if(!file.exists(outdir)) dir.create(outdir)

write.table(experiment_counts, paste0(outdir,'/countmat.txt'), quote=F, sep='\t')
write.table(data.frame(name=colnames(experiment_counts),sample=sprintf("sample_%02d",seq(1,ncol(experiment_counts))),stringsAsFactors=F), paste0(outdir,'/samples.txt'), quote=F, sep='\t', row.names=F)

simulate_experiment_countmat(smallseqfile, outdir=outdir,
                                  readmat=as.matrix(experiment_counts),
                                  fraglen=fraglen,
                                  fragsd=fragsd,
                                  readlen=readlen,
                                  error_rate=error_rate,
                                  paired=paired,
                                  seed=seed
)

unlink(smallseqfile)