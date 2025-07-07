#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

########################################
#TEMP WHILE WORKING ON SCRIPT
#args[1]<-"/PATH/REVAMP/outdir/blast_results" #working directory
#args[2]<-90 #blast query coverage cutoff
########################################
library(dplyr)
library(Biostrings)
setwd(as.character(args[1]))

all_results <- read.delim("ASV_blastn_nt.btab", header=FALSE, stringsAsFactors=FALSE)
colnames(all_results) <- c("ASV", "percent", "length", "taxid", "accession")

asv_seqs <- readDNAStringSet("../dada2/ASVs.fa")
asv_lengths <- data.frame(
  ASV = names(asv_seqs),
  asv_length = width(asv_seqs),
  stringsAsFactors = FALSE
)
all_results_joined <- left_join(all_results, asv_lengths, by = "ASV")

percent_cutoff <- as.numeric(args[2]) / 100
all_results_joined <- all_results_joined %>%
  mutate(min_required_length = asv_length * percent_cutoff)

all_results_trimmed <- all_results_joined %>%
  filter(length >= min_required_length)

all_results_trimmed$correction <- all_results_trimmed$percent / 100 * all_results_trimmed$length

all_results_best <- all_results_trimmed %>%
  group_by(ASV) %>%
  filter(percent == max(percent)) %>%
  summarise_all(~ toString(unique(.)))

all_results_best_original <- all_results_best %>%
  select("ASV", "percent", "length", "taxid", "accession", "correction")

write.table(all_results_best_original, file = 'ASV_blastn_nt_formatted.txt', sep = '\t', quote = FALSE, row.names = FALSE)
write.table(all_results_best, file = 'ASV_blastn_nt_formatted_verbose.txt', sep = '\t', quote = FALSE, row.names = FALSE)
