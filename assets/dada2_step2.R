#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

########################################
#TEMP WHILE WORKING ON SCRIPT
#args[1]<-"/Users/mcallister/Desktop/test_out/dada2" #working directory
#args[2]<-32000*0.7 #memory size Mb

########################################
library(dada2)

setwd(as.character(args[1]))

samples <- scan("../sample_order.txt", what="character")
sample_num <- length(samples)

filtered_forward_reads <- paste0(samples, "_R1_filtered.fq.gz")
filtered_reverse_reads <- paste0(samples, "_R2_filtered.fq.gz")

readstorun=450000*as.numeric(args[2])*0.7

message("DADA2 - Learning error Running")
err_forward_reads <- learnErrors(filtered_forward_reads, multithread = TRUE, nbases = readstorun)
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread = TRUE, nbases = readstorun)


fError <- plotErrors(err_forward_reads, nominalQ=TRUE)
pdf(file='errorFPlot.pdf', width=11, height = 8.5)
fError
dev.off()

rError <- plotErrors(err_reverse_reads, nominalQ=TRUE)
pdf(file='errorRPlot.pdf', width=11, height = 8.5)
rError
dev.off()

message("DADA2 - Dereplication Running")
derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
names(derep_forward) <- samples
derep_reverse <- derepFastq(filtered_reverse_reads, verbose=TRUE)
names(derep_reverse) <- samples

dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo", multithread = TRUE)
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pseudo", multithread = TRUE)

message("DADA2 - Merge Running")
merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse, derep_reverse, trimOverhang=TRUE, minOverlap=20)

message("DADA2 - Chimera Check Running")
seqtab <- makeSequenceTable(merged_amplicons)
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T)
message("Percent non-chimera:")
100*sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))

filtered_out <- read.delim("filtered_out_stats.txt", row.names=1, quote="")
summary_tab <- data.frame(row.names=samples, 
                          dada2_input=filtered_out[,1],
                          filtered=filtered_out[,2], 
                          dada_f=sapply(dada_forward, getN), 
                          dada_r=sapply(dada_reverse, getN), 
                          merged=sapply(merged_amplicons, getN), 
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
