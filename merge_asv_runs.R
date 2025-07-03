#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (!dir.exists(args[5])) {
  dir.create(args[5], recursive = TRUE)
}
args <- normalizePath(args, mustWork = TRUE)

#Note that this is designed to merge two sets of counts tables based on 
#identical sample names and merging identical ASVs. Unique samples and ASVs 
#will be appended to a new counts table. Requires ASVs_counts.tsv and ASVs.fa 
#from REVAMP. Output ASVs_counts.tsv, ASVs.fa, and ASV_nameChange.txt.

########################################
# TEMP WHILE WORKING ON SCRIPT
#args[1]<-"~/Desktop/OSU_sequencing_analysis/temp_REVAMP/20230331_MiSeq1040_Ax20_Pool2_Lane1/16S_out/dada2/ASVs_counts.tsv" #path to ASV count table 1
#args[2]<-"~/Desktop/OSU_sequencing_analysis/temp_REVAMP/20230331_MiSeq1040_Ax20_Pool2_Lane1/16S_out/dada2/ASVs.fa" #path to ASV fasta file 1
#args[3]<-"~/Desktop/OSU_sequencing_analysis/temp_REVAMP/20230905_MiSeq1079_Ax20_Pool2_Lane2/AX_16S_out/dada2/ASVs_counts.tsv" #path to ASV count table 2
#args[4]<-"~/Desktop/OSU_sequencing_analysis/temp_REVAMP/20230905_MiSeq1079_Ax20_Pool2_Lane2/AX_16S_out/dada2/ASVs.fa" #path to ASV fasta file 2
#args[5]<-"~/Desktop/OSU_sequencing_analysis/temp_REVAMP/combine_Ax_lanes_16S" #outdir
########################################
library("dplyr")
library("tidyr")
library("stringr")
library("Biostrings")

setwd(args[5])

#Read in the Pair 1 Counts/Sequences
asv_count1 <- read.delim(as.character(args[1]), header=TRUE, stringsAsFactors=FALSE)
names(asv_count1)[names(asv_count1) == "x"] <- "ASVs"

fasta_file1 <- readDNAStringSet(as.character(args[2]))
fasta1_df <- data.frame(
  ASVs = names(fasta_file1),
  sequence = as.character(fasta_file1),
  stringsAsFactors = FALSE
)

merge1 <- asv_count1 %>% left_join(fasta1_df, by = "ASVs")

#Read in the Pair 2 Counts/Sequences
asv_count2 <- read.delim(as.character(args[3]), header=TRUE, stringsAsFactors=FALSE)
names(asv_count2)[names(asv_count2) == "x"] <- "ASVs"

fasta_file2 <- readDNAStringSet(as.character(args[4]))
fasta2_df <- data.frame(
  ASVs = names(fasta_file2),
  sequence = as.character(fasta_file2),
  stringsAsFactors = FALSE
)

merge2 <- asv_count2 %>% left_join(fasta2_df, by = "ASVs")

#Merge the two count tables
merged <- merge1 %>%
  full_join(merge2, by = "sequence", suffix = c(".merge1", ".merge2"))

names(merged)[names(merged) == "ASVs.merge1"] <- "ASVs_merge1"
names(merged)[names(merged) == "ASVs.merge2"] <- "ASVs_merge2"

#Combine counts for overlapping sample columns
sample_cols_merge1 <- grep("\\.merge1$", names(merged), value = TRUE)
sample_cols_merge2 <- grep("\\.merge2$", names(merged), value = TRUE)
sample_cols <- unique(sub("\\.merge[12]$", "", c(sample_cols_merge1, sample_cols_merge2)))

for (col in sample_cols) {
  merged[[col]] <- rowSums(
    cbind(
      merged[[paste0(col, ".merge1")]],
      merged[[paste0(col, ".merge2")]]
    ),
    na.rm = TRUE
  )
}

# Save original ASV names (left = dataset 1; right = dataset 2)
merged <- merged %>%
  mutate(
    ASVs_original = paste0(
      coalesce(ASVs_merge1, "NA"), ";",
      coalesce(ASVs_merge2, "NA")
    )
  ) %>%
  select(sequence, ASVs_original, starts_with("MP_")) %>%
  select(-ends_with(".merge1"), -ends_with(".merge2"))

#Compute total counts for ordering and renumber ASVs
merged <- merged %>%
  mutate(across(starts_with("MP_"), ~replace_na(., 0)))

merged <- merged %>%
  mutate(total_counts = rowSums(across(starts_with("MP_")))) %>%
  arrange(desc(total_counts)) %>%
  mutate(ASVs_new = paste0("ASV_", row_number())) %>%
  select(ASVs_new, sequence, ASVs_original, everything(), -total_counts)

#Separate to outputs
ASV_counts_df <- merged %>% select(-sequence, -ASVs_original)
names(ASV_counts_df)[names(ASV_counts_df) == "ASVs_new"] <- "x"
write.table(ASV_counts_df, 
            file = "ASVs_counts.tsv", 
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

asvs <- merged$ASVs_new
sequences <- merged$sequence
fasta_sequences <- DNAStringSet(sequences)
names(fasta_sequences) <- asvs
writeXStringSet(fasta_sequences, filepath = "ASVs.fa")

ASV_name_change <- merged %>% select(ASVs_new, ASVs_original)
write.table(ASV_name_change, 
            file = "ASVs_nameChange.txt", 
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
