#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

########################################
#TEMP WHILE WORKING ON SCRIPT
# args[1]<-"/Users/mcallister/Desktop/sean_test/dada2" #working directory
# args[2]<-32000*0.7 #memory size Mb
# args[3]<-"forward" #which read direction to keep

########################################
library(dada2)

setwd(as.character(args[1]))

samples <- scan("../sample_order.txt", what="character")
sample_num <- length(samples)

readstorun=450000*as.numeric(args[2])*0.7

if (as.character(args[3]) == "forward") {
  filtered_forward_reads <- paste0(samples, "_R1_filtered.fq.gz")

  message("DADA2 - Learning error Running")
  err_forward_reads <- learnErrors(filtered_forward_reads, multithread = TRUE, nbases = readstorun)
  
  message("DADA2 - Dereplication Running")
  derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
  names(derep_forward) <- samples
  
  dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo", multithread = TRUE)
  
  message("DADA2 - Chimera Check Running")
  seqtab <- makeSequenceTable(dada_forward)
  seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T)
  message("Percent non-chimera:")
  100*sum(seqtab.nochim)/sum(seqtab)
  
  getN <- function(x) sum(getUniques(x))
  
  filtered_out <- read.delim("filtered_out_stats.txt", row.names=1, quote="")
  
  summary_tab <- data.frame(row.names=samples, 
                            dada2_input=filtered_out[,1],
                            filtered=filtered_out[,2], 
                            dada_f=sapply(dada_forward, getN), 
                            nonchim=rowSums(seqtab.nochim), 
                            final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))
  write.table(summary_tab, file='ReadTrimSummary.txt', sep= '\t', quote=FALSE)
  
  message("DADA2 - Creating ASV outs")
  asv_seqs <- colnames(seqtab.nochim)
  asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
  for (i in 1:dim(seqtab.nochim)[2]) {
    asv_headers[i] <- paste(">ASV", i, sep="_")
  }
  asv_fasta <- c(rbind(asv_headers, asv_seqs))
  write(asv_fasta, "ASVs.fa")
  
  asv_tab <- t(seqtab.nochim)
  row.names(asv_tab) <- sub(">", "", asv_headers)
  write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)
}


if (as.character(args[3]) == "reverse") {
  filtered_reverse_reads <- paste0(samples, "_R2_filtered.fq.gz")
  
  message("DADA2 - Learning error Running")
  err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread = TRUE, nbases = readstorun)
  
  message("DADA2 - Dereplication Running")
  derep_reverse <- derepFastq(filtered_reverse_reads, verbose=TRUE)
  names(derep_reverse) <- samples
  
  dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pseudo", multithread = TRUE)
  
  message("DADA2 - Chimera Check Running")
  seqtab <- makeSequenceTable(dada_reverse)
  seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T)
  message("Percent non-chimera:")
  100*sum(seqtab.nochim)/sum(seqtab)
  
  getN <- function(x) sum(getUniques(x))
  
  filtered_out <- read.delim("filtered_out_stats.txt", row.names=1, quote="")
  summary_tab <- data.frame(row.names=samples, 
                            dada2_input=filtered_out[,1],
                            filtered=filtered_out[,2], 
                            dada_r=sapply(dada_reverse, getN), 
                            nonchim=rowSums(seqtab.nochim), 
                            final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))
  write.table(summary_tab, file='ReadTrimSummary.txt', sep= '\t', quote=FALSE)
  
  message("DADA2 - Creating ASV outs")
  asv_seqs <- colnames(seqtab.nochim)
  asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
  for (i in 1:dim(seqtab.nochim)[2]) {
    asv_headers[i] <- paste(">ASV", i, sep="_")
  }
  asv_fasta <- c(rbind(asv_headers, asv_seqs))
  write(asv_fasta, "ASVs.fa")
  
  asv_tab <- t(seqtab.nochim)
  row.names(asv_tab) <- sub(">", "", asv_headers)
  write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)
}
