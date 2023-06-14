#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

########################################
#TEMP WHILE WORKING ON SCRIPT
#args[1]<-"/PATH/REVAMP/outdir/dada2" #working directory
#args[2]<-120 #minlength 
#args[3]<-TRUE #rm.phix
#args[4]<-2 #truncQ
#args[5]<-2 #maxEE No .1
#args[6]<-2 #maxEE No. 2
#args[7]<-0 #trimRight
#args[8]<-0 #trimLeft
########################################
library(dada2)

setwd(as.character(args[1]))

samples <- scan("../sample_order.txt", what="character")
sample_num <- length(samples)

forward_reads <- paste0(samples, "_R1_trimmed.fq.gz")
reverse_reads <- paste0(samples, "_R2_trimmed.fq.gz")

filtered_forward_reads <- paste0(samples, "_R1_filtered.fq.gz")
filtered_reverse_reads <- paste0(samples, "_R2_filtered.fq.gz")

setwd(paste0(as.character(args[1]), "/../cutadapt"))
rawFQuality <- plotQualityProfile(forward_reads)
setwd(as.character(args[1]))
pdf(file='rawFQualityPlot.pdf', width=11, height = 7*(sample_num/7))
print(rawFQuality)
dev.off()

setwd(paste0(as.character(args[1]), "/../cutadapt"))
rawRQuality <- plotQualityProfile(reverse_reads)
setwd(as.character(args[1]))
pdf(file='rawRQualityPlot.pdf', width=11, height = 7*(sample_num/7))
print(rawRQuality)
dev.off()

setwd(paste0(as.character(args[1]), "/../cutadapt"))
filtered_out <- filterAndTrim(fwd = forward_reads, 
                              filt = paste0(args[1], "/", filtered_forward_reads), 
                              rev = reverse_reads, 
                              filt.rev = paste0(args[1], "/", filtered_reverse_reads),
                              multithread = TRUE,
                              minLen = as.numeric(args[2]),
                              rm.phix= as.logical(args[3]),
                              truncQ = as.numeric(args[4]), 
                              maxEE = c(as.numeric(args[5]), as.numeric(args[6])),
                              trimRight = as.numeric(args[7]),
                              trimLeft = as.numeric(args[8]))

setwd(as.character(args[1]))
write.table(filtered_out, file='filtered_out_stats.txt', sep= '\t', quote=FALSE)

trimFQuality <- plotQualityProfile(filtered_forward_reads)
pdf(file='trimFQualityPlot.pdf', width=11, height = 7*(sample_num/7))
print(trimFQuality)
dev.off()

trimRQuality <- plotQualityProfile(filtered_reverse_reads)
pdf(file='trimRQualityPlot.pdf', width=11, height = 7*(sample_num/7))
print(trimRQuality)
dev.off()
