#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

seed <- as.integer(args[1])
if(is.na(seed)) seed <- NULL

set.seed(seed)

fams <- c("HERVK", "HERVW", "HERVH", "HERVE")
candidates = read.table('candidates.txt', sep='\t', header=T, stringsAsFactors=F)

locnames <- data.frame(do.call(cbind,lapply(fams,function(f){sample(candidates$locus[candidates$family==f],100,replace=T)})),stringsAsFactors=F)
proportions <- data.frame(do.call(cbind,lapply(1:100,function(i){ifelse(candidates$locus %in% locnames[i,],1,0)})))
row.names(proportions) <- candidates$locus
names(proportions) <- sprintf('sample_%02d',seq(1,100,1))
proportions <- proportions / colSums(proportions)
write.table(proportions,'simulation_proportions.txt', sep='\t', quote=F)
