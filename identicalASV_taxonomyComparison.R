#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

########################################
#TEMP WHILE WORKING ON SCRIPT
#args[1]<-"/PATH/to/REVAMP/outdir/1" #working directory comp 1
#args[2]<-"/PATH/to/REVAMP/outdir/2" #working directory comp 2
#args[3]<-"/PATH/to/comparison/outdir" #path to outdir
########################################
library(dplyr)
library(Biostrings)

setwd(as.character(args[1]))
taxonomy_file <- list.files(path = "./ASV2Taxonomy", pattern = "_asvTaxonomyTable\\.txt$", full.names = TRUE)
asv_tax_strings_1 <- read.delim(as.character(taxonomy_file), header=TRUE, stringsAsFactors=FALSE)
asv_tax_strings_1 <- asv_tax_strings_1 %>%
  mutate(taxonomy_str = apply(select(., -ASV), 1, function(x) paste(x, collapse = ";")))
asv_ids_1 <- read.delim("./blast_results/ASV_blastn_nt_formatted.txt", header=TRUE, stringsAsFactors=FALSE)
asv_seqs_1 <- readDNAStringSet("./dada2/ASVs.fa")
asv_lengths_1 <- data.frame(
  ASV = names(asv_seqs_1),
  asv_length = width(asv_seqs_1),
  stringsAsFactors = FALSE
)

setwd(as.character(args[2]))
taxonomy_file <- list.files(path = "./ASV2Taxonomy", pattern = "_asvTaxonomyTable\\.txt$", full.names = TRUE)
asv_tax_strings_2 <- read.delim(as.character(taxonomy_file), header=TRUE, stringsAsFactors=FALSE)
asv_tax_strings_2 <- asv_tax_strings_2 %>%
  mutate(taxonomy_str = apply(select(., -ASV), 1, function(x) paste(x, collapse = ";")))
asv_ids_2 <- read.delim("./blast_results/ASV_blastn_nt_formatted.txt", header=TRUE, stringsAsFactors=FALSE)
asv_seqs_2 <- readDNAStringSet("./dada2/ASVs.fa")
asv_lengths_2 <- data.frame(
  ASV = names(asv_seqs_2),
  asv_length = width(asv_seqs_2),
  stringsAsFactors = FALSE
)

asv_seqs_1_sorted <- asv_seqs_1[order(names(asv_seqs_1))]
asv_seqs_2_sorted <- asv_seqs_2[order(names(asv_seqs_2))]
if (!identical(names(asv_seqs_1_sorted), names(asv_seqs_2_sorted)) || 
    !identical(as.character(asv_seqs_1_sorted), as.character(asv_seqs_2_sorted))) {
  stop("ERROR: Only identical ASVs between runs can be compared.")
}

asv_comparison_results <- data.frame(ASV = names(asv_seqs_1), stringsAsFactors = FALSE) %>%
  left_join(asv_lengths_1, by = "ASV") %>%
  left_join(asv_tax_strings_1 %>% select(ASV, taxa_string_1 = taxonomy_str), by = "ASV") %>%
  left_join(asv_ids_1 %>% select(ASV, match_pid_1 = percent, match_length_1 = length), by = "ASV") %>%
  left_join(asv_tax_strings_2 %>% select(ASV, taxa_string_2 = taxonomy_str), by = "ASV") %>%
  left_join(asv_ids_2 %>% select(ASV, match_pid_2 = percent, match_length_2 = length), by = "ASV")

get_match_type <- function(taxa1, taxa2) {
  tax1 <- unlist(strsplit(taxa1, ";"))
  tax2 <- unlist(strsplit(taxa2, ";"))
  
  tax1 <- rep(tax1, length.out = 7)
  tax2 <- rep(tax2, length.out = 7)
  
  unknown_strings <- c("Unknown;NA;NA;NA;NA;NA;NA", "Environmental Unknown;NA;NA;NA;NA;NA;NA")
  
  if (taxa1 %in% unknown_strings && !(taxa2 %in% unknown_strings)) {
    return("recovered_from_unknown_1to2")
  } else if (taxa2 %in% unknown_strings && !(taxa1 %in% unknown_strings)) {
    return("recovered_from_unknown_2to1")
  }
  
  if (identical(tax1, tax2)) return("identical")
  
  levels <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  for (i in seq_along(tax1)) {
    if (!identical(tax1[i], tax2[i])) {
      if (i == 1) return("no_match")
      return(paste0("same_", levels[i - 1]))
    }
  }
  
  return("unknown_match")
}

asv_comparison_results <- asv_comparison_results %>%
  rowwise() %>%
  mutate(match_type = get_match_type(taxa_string_1, taxa_string_2)) %>%
  ungroup()

setwd(as.character(args[3]))
write.table(asv_comparison_results, file = 'comparison_ASV_results.txt', sep = '\t', quote = FALSE, row.names = FALSE)

total_asvs <- nrow(asv_comparison_results)
match_counts <- table(asv_comparison_results$match_type)
format_pct <- function(type) {
  count <- if (type %in% names(match_counts)) match_counts[[type]] else 0
  pct <- 100 * count / total_asvs
  if (pct < 0.01 && pct > 0) {
    return("<0.01%")
  } else {
    return(sprintf("%.2f%%", pct))
  }
}
cat("Comparison Results:\n")
cat("# ASVs:", total_asvs, "\n")
cat("% identical matches:", format_pct("identical"), "\n")
cat("% recovered unknown 1to2 matches:", format_pct("recovered_from_unknown_1to2"), "\n")
cat("% recovered unknown 2to1 matches:", format_pct("recovered_from_unknown_2to1"), "\n")
cat("% same Genus matches:", format_pct("same_genus"), "\n")
cat("% same Family matches:", format_pct("same_family"), "\n")
higher_levels <- c("same_order", "same_class", "same_phylum", "same_kingdom")
sum_pct <- 100 * sum(match_counts[higher_levels], na.rm = TRUE) / total_asvs
if (sum_pct < 0.01 && sum_pct > 0) {
  higher_result <- "<0.01%"
} else {
  higher_result <- sprintf("%.2f%%", sum_pct)
}
cat("% same higher taxonomy matches:", higher_result, "\n")
