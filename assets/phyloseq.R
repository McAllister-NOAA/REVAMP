#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

########################################
#TEMP WHILE WORKING ON SCRIPT
# args[1]<-"/Users/mcallister/Desktop/Angie/COI_27Aug21/MTeDNA/MTeDNA_out/Figures" #FIGURE OUT directory
# args[2]<-"/Users/mcallister/Desktop/Angie/COI_27Aug21/MTeDNA" #Working directory
# args[3]<-"MTeDNA_out" #outdirectory name
# args[4]<-"TRUE" #control flag
# args[5]<-"FALSE" #filter low quality samples flag
# args[6]<-"FALSE" #replicateFlag
# args[7]<-"TRUE" #sitelabelFlag
# args[8]<-"5" #percent filter for taxa
# args[9]<-"FALSE" #FILTER NAs where appropriate
# args[10]<-"FALSE" #whether taxa of interest file is given
# args[11]<-"TRUE" #groups defined in sample metadata file
# args[12]<-"3" #number of groups defined in sample metadata file
# args[13]<-"Order" #category to filter on for taxa of interest
# args[14]<-"/Users/mcallister/Desktop/chris_test/18S/choice_taxa.txt" #location of taxa of interest one per line
# args[15]<-"FALSE" #whether chem data has been provided
# args[16]<-"NULL" #location of chem header file (one per line)
########################################
library("ggplot2")
library("dplyr")
library("phyloseq")
#library("phylosmith")
#library("data.table")
library("ggpubr")
library("ggalt")
library("ggrepel")
#library("plyr")
library("vegan")
library("spatstat")

setwd(as.character(args[1]))
theme_set(theme_bw())

controlFlag <- as.logical(args[4])
filteredLowQualSamples <- as.logical(args[5])
replicateFlag <- as.logical(args[6])
sitelabelFlag <- as.logical(args[7])
NAremoveFlag <- as.logical(args[9])
taxaofinterestFlag <- as.logical(args[10])
groupingFlag <- as.logical(args[11])
numberofGroups <- as.numeric(args[12])
chemDataFlag <- as.logical(args[15])

filter_percent <- as.numeric(args[8]) #Taxa below this were filtered to zzOther

if (chemDataFlag == TRUE) {
chem_headers <- read.delim(as.character(args[16]), header=FALSE, stringsAsFactors=FALSE)
chem_headers <- as.vector(chem_headers)
}

##################################
#
#  Import ASV-based, raw reads, no qual filter (phylo1)
#
##################################
asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/ASV2Taxonomy/ASVs_counts_NOUNKNOWNS.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
row.names(asv_count) <- asv_count$x
asv_count <- asv_count %>% select(-x)
asv_count_mat <- as.matrix(asv_count)

asv_taxonomy <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/ASV2Taxonomy/",as.character(args[3]),"_asvTaxonomyTable_NOUNKNOWNS.txt", sep = ""), header=TRUE, stringsAsFactors=FALSE)
row.names(asv_taxonomy) <- asv_taxonomy$ASV
asv_taxonomy <- asv_taxonomy %>% select(-ASV)
asv_taxonomy_mat <- as.matrix(asv_taxonomy)

sample_metadata <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/sample_metadata_forR.txt", sep = ""), header=TRUE, stringsAsFactors=TRUE)
row.names(sample_metadata) <- sample_metadata$Sample
sample_metadata <- sample_metadata %>% select(-Sample)

ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
TAX = tax_table(asv_taxonomy_mat)
samples = sample_data(sample_metadata)

phylo1 <- phyloseq(ASV, TAX, samples)

##################################
#
#  Import ASV-based, raw reads, qual filtered (phylo2)
#
##################################
if (controlFlag == TRUE && filteredLowQualSamples == TRUE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_lowEffortSamplesRemoved.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo2 <- phyloseq(ASV, TAX, samples)
} else if (controlFlag == TRUE && filteredLowQualSamples == FALSE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_controlsRemoved.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo2 <- phyloseq(ASV, TAX, samples)
} else if (controlFlag == FALSE && filteredLowQualSamples == TRUE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_lowEffortSamplesRemoved.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo2 <- phyloseq(ASV, TAX, samples)
} else if (controlFlag == FALSE && filteredLowQualSamples == FALSE) {
  phylo2 <- phylo1
}

##################################
#
#  Import ASV-based, relative abund, qual filtered (phylo3)
#
##################################
asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_percentabund.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
row.names(asv_count) <- asv_count$x
asv_count <- asv_count %>% select(-x)
asv_count_mat <- as.matrix(asv_count)
ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
phylo3 <- phyloseq(ASV, TAX, samples)


##################################
#
#  Import Taxonomy-based, raw reads, qual filtered (phylo10)
#
##################################
if (controlFlag == TRUE && filteredLowQualSamples == TRUE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_lowEffortSamplesRemoved.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo10 <- phyloseq(ASV, TAX, samples)
} else if (controlFlag == TRUE && filteredLowQualSamples == FALSE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_controlsRemoved.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo10 <- phyloseq(ASV, TAX, samples)
} else if (controlFlag == FALSE && filteredLowQualSamples == TRUE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_lowEffortSamplesRemoved.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo10 <- phyloseq(ASV, TAX, samples)
} else if (controlFlag == FALSE && filteredLowQualSamples == FALSE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy_NOUNKNOWNS.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo10 <- phyloseq(ASV, TAX, samples)
}

##################################
#
#  Import Taxonomy-based, relative abund, qual filtered (phylo11)
#
##################################
asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
row.names(asv_count) <- asv_count$x
asv_count <- asv_count %>% select(-x)
asv_count_mat <- as.matrix(asv_count)
ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
phylo11 <- phyloseq(ASV, TAX, samples)

##################################
#
#  Import ASV-based, relative abund, qual filtered, grouped on replicates/sites (phylo7a/b)
#
##################################
#ADD GROUPED SAMPLE METADATA
if (replicateFlag == TRUE) {
  grouped_sample_metadata <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/sample_metadata_NOUNKNOWNS_percentabund_groupedByReplicates.tsv", sep = ""), header=TRUE, stringsAsFactors=TRUE)
  grouped_sample_metadata <- grouped_sample_metadata %>% filter(!is.na(Sample))
  row.names(grouped_sample_metadata) <- grouped_sample_metadata$Sample
  grouped_sample_metadata <- grouped_sample_metadata %>% select(-Sample)
  grouped_samplesA = sample_data(grouped_sample_metadata)
  
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_percentabund_groupedByReplicates.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo7a <- phyloseq(ASV, TAX, grouped_samplesA)
}
if (sitelabelFlag == TRUE) {
  grouped_sample_metadata <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/sample_metadata_NOUNKNOWNS_percentabund_groupedBySites.tsv", sep = ""), header=TRUE, stringsAsFactors=TRUE)
  grouped_sample_metadata <- grouped_sample_metadata %>% filter(!is.na(Sample))
  row.names(grouped_sample_metadata) <- grouped_sample_metadata$Sample
  grouped_sample_metadata <- grouped_sample_metadata %>% select(-Sample)
  grouped_samplesB = sample_data(grouped_sample_metadata)
  
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_percentabund_groupedBySites.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo7b <- phyloseq(ASV, TAX, grouped_samplesB)
}

##################################
#
#  Import Taxonomy-based, relative abund, qual filtered, grouped on replicates/sites (phylo15a/b)
#
##################################
if (replicateFlag == TRUE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund_groupedByReplicates.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo15a <- phyloseq(ASV, TAX, grouped_samplesA)
}
if (sitelabelFlag == TRUE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund_groupedBySites.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo15b <- phyloseq(ASV, TAX, grouped_samplesB)
}

##################################
#
#  Import ASV-based, raw reads, qual filtered, taxa filtered (phylo4)
#
##################################
#CHANGE TO NEW FILTERED TAXONOMY FILE
asv_taxonomy <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt", sep = ""), header=TRUE, stringsAsFactors=FALSE)
row.names(asv_taxonomy) <- asv_taxonomy$ASV
asv_taxonomy <- asv_taxonomy %>% select(-ASV)
asv_taxonomy_mat <- as.matrix(asv_taxonomy)
TAX = tax_table(asv_taxonomy_mat)

if (controlFlag == TRUE && filteredLowQualSamples == TRUE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_lowEffortSamplesRemoved.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo4 <- phyloseq(ASV, TAX, samples)
} else if (controlFlag == TRUE && filteredLowQualSamples == FALSE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_controlsRemoved.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo4 <- phyloseq(ASV, TAX, samples)
} else if (controlFlag == FALSE && filteredLowQualSamples == TRUE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_lowEffortSamplesRemoved.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo4 <- phyloseq(ASV, TAX, samples)
} else if (controlFlag == FALSE && filteredLowQualSamples == FALSE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/ASV2Taxonomy/ASVs_counts_NOUNKNOWNS.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo4 <- phyloseq(ASV, TAX, samples)
}

##################################
#
#  Import ASV-based, relative abund, qual filtered, taxa filtered (phylo5)
#
##################################
asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_percentabund.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
row.names(asv_count) <- asv_count$x
asv_count <- asv_count %>% select(-x)
asv_count_mat <- as.matrix(asv_count)
ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
phylo5 <- phyloseq(ASV, TAX, samples)

##################################
#
#  Import Taxonomy-based, raw reads, qual filtered, taxa filtered (phylo12)
#
##################################
if (controlFlag == TRUE && filteredLowQualSamples == TRUE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_lowEffortSamplesRemoved.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo12 <- phyloseq(ASV, TAX, samples)
} else if (controlFlag == TRUE && filteredLowQualSamples == FALSE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_controlsRemoved.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo12 <- phyloseq(ASV, TAX, samples)
} else if (controlFlag == FALSE && filteredLowQualSamples == TRUE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_lowEffortSamplesRemoved.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo12 <- phyloseq(ASV, TAX, samples)
} else if (controlFlag == FALSE && filteredLowQualSamples == FALSE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy_NOUNKNOWNS.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo12 <- phyloseq(ASV, TAX, samples)
}

##################################
#
#  Import Taxonomy-based, relative abund, qual filtered, taxa filtered (phylo13)
#
##################################
asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
row.names(asv_count) <- asv_count$x
asv_count <- asv_count %>% select(-x)
asv_count_mat <- as.matrix(asv_count)
ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
phylo13 <- phyloseq(ASV, TAX, samples)

##################################
#
#  Import ASV-based, relative abund, qual filtered, taxa filtered, grouped on replicates/sites (phylo9a/b)
#
##################################
if (replicateFlag == TRUE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_percentabund_groupedByReplicates.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo9a <- phyloseq(ASV, TAX, grouped_samplesA)
}
if (sitelabelFlag == TRUE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_percentabund_groupedBySites.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo9b <- phyloseq(ASV, TAX, grouped_samplesB)
}

##################################
#
#  Import Taxonomy-based, relative abund, qual filtered, taxa filtered, grouped on replicates/sites (phylo17a/b)
#
##################################
if (replicateFlag == TRUE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund_groupedByReplicates.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo17a <- phyloseq(ASV, TAX, grouped_samplesA)
}
if (sitelabelFlag == TRUE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund_groupedBySites.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo17b <- phyloseq(ASV, TAX, grouped_samplesB)
}


message("#######################################################################################")
message("#")
message("# Phyloseq objects have been created to work with:")
message("#")
message("# ASV-based, raw reads, no qual filter (phylo1) ALL")
message("# ASV-based, raw reads, qual filtered (phylo2) ALL")
message("# ASV-based, relative abund, qual filtered (phylo3) ALL")
message("# ASV-based, raw reads, qual filtered, taxa filtered (phylo4) ALL")
message("# ASV-based, relative abund, qual filtered, taxa filtered (phylo5) ALL")
message("# phylo6 depracated")
message("# ASV-based, relative abund, qual filtered, grouped on replicates/sites (phylo7a/b) IF replicates/sites called")
message("# phylo8 depracated")
message("# ASV-based, relative abund, qual filtered, taxa filtered, grouped on replicates/sites (phylo9a/b) IF replicates/sites called")
message("#")
message("# Taxonomy-based, raw reads, qual filtered (phylo10) ALL")
message("# Taxonomy-based, relative abund, qual filtered (phylo11) ALL")
message("# Taxonomy-based, raw reads, qual filtered, taxa filtered (phylo12) ALL")
message("# Taxonomy-based, relative abund, qual filtered, taxa filtered (phylo13) ALL")
message("# phylo14 depracated")
message("# Taxonomy-based, relative abund, qual filtered, grouped on replicates/sites (phylo15a/b) IF replicates/sites called")
message("# phylo16 depracated")
message("# Taxonomy-based, relative abund, qual filtered, taxa filtered, grouped on replicates/sites (phylo17a/b) IF replicates/sites called")
message("#")
message("#######################################################################################")

sample_order <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/sample_order.txt", sep = ""), header=FALSE, stringsAsFactors=FALSE)
sample_order <- as.vector(sample_order)

if (replicateFlag == TRUE) {
  replicateorder <- as.data.frame(sample_metadata$replicates)
  replicateorder <- distinct(replicateorder)
  replicateorder <- as.vector(replicateorder)
}
if (sitelabelFlag == TRUE) {
  siteorder <- as.data.frame(sample_metadata$sites)
  siteorder <- distinct(siteorder)
  siteorder <- as.vector(siteorder)
}

########################
#
# BARCHARTS
#
########################
tryCatch({
setwd(paste0(as.character(args[1]), "/02_Barcharts/read_count", sep = ""))

base_Phylum_count <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/ASV2Taxonomy/",as.character(args[3]),"_asvTaxonomyTable_NOUNKNOWNS.txt | cut -f3 | grep -v ","\"","Phylum","\""," | sort | uniq | wc -l", sep = ""),
                                  intern = TRUE))
base_Class_count <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/ASV2Taxonomy/",as.character(args[3]),"_asvTaxonomyTable_NOUNKNOWNS.txt | cut -f4 | grep -v ","\"","Class","\""," | sort | uniq | wc -l", sep = ""),
                                 intern = TRUE))
base_Order_count <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/ASV2Taxonomy/",as.character(args[3]),"_asvTaxonomyTable_NOUNKNOWNS.txt | cut -f5 | grep -v ","\"","Order","\""," | sort | uniq | wc -l", sep = ""),
                                 intern = TRUE))
base_Family_count <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/ASV2Taxonomy/",as.character(args[3]),"_asvTaxonomyTable_NOUNKNOWNS.txt | cut -f6 | grep -v ","\"","Family","\""," | sort | uniq | wc -l", sep = ""),
                                  intern = TRUE))
base_Genus_count <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/ASV2Taxonomy/",as.character(args[3]),"_asvTaxonomyTable_NOUNKNOWNS.txt | cut -f7 | grep -v ","\"","Genus","\""," | sort | uniq | wc -l", sep = ""),
                                 intern = TRUE))
base_Species_count <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/ASV2Taxonomy/",as.character(args[3]),"_asvTaxonomyTable_NOUNKNOWNS.txt | cut -f8 | grep -v ","\"","Species","\""," | sort | uniq | wc -l", sep = ""),
                                   intern = TRUE))

tryCatch({
for (i in c("Phylum")) {
  if (NAremoveFlag == TRUE) {
    bar_plot <- plot_bar(subset_taxa(phylo1, eval(parse(text = paste0("!is.na(",i,")")))), fill = i) +
      geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
      ylab("Abundance (read count)") +
      theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
  } else {
    bar_plot <- plot_bar(phylo1, fill = i) +
      geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
      ylab("Abundance (read count)") +
      theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
  }
  bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = sample_order$V1)
  bar_plot_legend <- get_legend(bar_plot)
  bar_plot <- bar_plot + theme(legend.position='none')
  legend_call <- eval(parse(text = paste0("base_",i,"_count")))
  w_choice <- 6.4 + ((legend_call*2.3)/20)
  pdf(file=paste0('barplot_readcount_allsamples_alltaxa_',i,'_legend.pdf', sep = ""), width = w_choice, height = 6)
  print(as_ggplot(bar_plot_legend))
  dev.off()
  pdf(file=paste0('barplot_readcount_allsamples_alltaxa_',i,'.pdf', sep = ""), width = 11, height = 8.5)
  print(bar_plot)
  dev.off()
}
}, error=function(e){})

tryCatch({
if (controlFlag == TRUE || filteredLowQualSamples == TRUE) {
  for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    if (NAremoveFlag == TRUE) {
      bar_plot <- plot_bar(subset_taxa(phylo2, eval(parse(text = paste0("!is.na(",i,")")))), fill = i) +
        geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
        ylab("Abundance (read count)") +
        theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
    } else {
      bar_plot <- plot_bar(phylo2, fill = i) +
        geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
        ylab("Abundance (read count)") +
        theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
    }
    bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = sample_order$V1)
    bar_plot_legend <- get_legend(bar_plot)
    bar_plot <- bar_plot + theme(legend.position='none')
    legend_call <- eval(parse(text = paste0("base_",i,"_count")))
    w_choice <- 6.4 + ((legend_call*2.3)/20)
    pdf(file=paste0('barplot_readcount_filtsamples_alltaxa_',i,'_legend.pdf', sep = ""), width = w_choice, height = 6)
    print(as_ggplot(bar_plot_legend))
    dev.off()
    pdf(file=paste0('barplot_readcount_filtsamples_alltaxa_',i,'.pdf', sep = ""), width = 11, height = 8.5)
    print(bar_plot)
    dev.off()
  }
} else {
  for (i in c("Class", "Order", "Family", "Genus", "Species")) {
    if (NAremoveFlag == TRUE) {
      bar_plot <- plot_bar(subset_taxa(phylo1, eval(parse(text = paste0("!is.na(",i,")")))), fill = i) +
        geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
        ylab("Abundance (read count)") +
        theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
    } else {
      bar_plot <- plot_bar(phylo1, fill = i) +
        geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
        ylab("Abundance (read count)") +
        theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
    }
    bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = sample_order$V1)
    bar_plot_legend <- get_legend(bar_plot)
    bar_plot <- bar_plot + theme(legend.position='none')
    legend_call <- eval(parse(text = paste0("base_",i,"_count")))
    w_choice <- 6.4 + ((legend_call*2.3)/20)
    pdf(file=paste0('barplot_readcount_allsamples_alltaxa_',i,'_legend.pdf', sep = ""), width = w_choice, height = 6)
    print(as_ggplot(bar_plot_legend))
    dev.off()
    pdf(file=paste0('barplot_readcount_allsamples_alltaxa_',i,'.pdf', sep = ""), width = 11, height = 8.5)
    print(bar_plot)
    dev.off()
  }
}
}, error=function(e){})

#Filtered taxa barcharts
tryCatch({
Phylum_count <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f3 | grep -v ","\"","Phylum","\""," | sort | uniq | wc -l", sep = ""),
                                  intern = TRUE))
Phylum_NA <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f3 | grep ","\"","^NA$","\""," | sort | uniq | wc -l", sep = ""),
                               intern = TRUE))
Class_count <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f4 | grep -v ","\"","Class","\""," | sort | uniq | wc -l", sep = ""),
                                 intern = TRUE))
Class_NA <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f4 | grep ","\"","^NA$","\""," | sort | uniq | wc -l", sep = ""),
                                 intern = TRUE))
Order_count <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f5 | grep -v ","\"","Order","\""," | sort | uniq | wc -l", sep = ""),
                                 intern = TRUE))
Order_NA <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f5 | grep ","\"","^NA$","\""," | sort | uniq | wc -l", sep = ""),
                                 intern = TRUE))
Family_count <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f6 | grep -v ","\"","Family","\""," | sort | uniq | wc -l", sep = ""),
                                  intern = TRUE))
Family_NA <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f6 | grep ","\"","^NA$","\""," | sort | uniq | wc -l", sep = ""),
                                  intern = TRUE))
Genus_count <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f7 | grep -v ","\"","Genus","\""," | sort | uniq | wc -l", sep = ""),
                                 intern = TRUE))
Genus_NA <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f7 | grep ","\"","^NA$","\""," | sort | uniq | wc -l", sep = ""),
                                 intern = TRUE))
Species_count <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f8 | grep -v ","\"","Species","\""," | sort | uniq | wc -l", sep = ""),
                                   intern = TRUE))
Species_NA <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f8 | grep ","\"","^NA$","\""," | sort | uniq | wc -l", sep = ""),
                                   intern = TRUE))
gg_color <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
  count_call <- eval(parse(text = paste0(i,"_count")))
  NA_call <- eval(parse(text = paste0(i,"_NA")))
  temp_col_pallete <- c(gg_color(count_call-NA_call-1), "grey")
  if (NAremoveFlag == TRUE) {
    bar_plot <- plot_bar(subset_taxa(phylo4, eval(parse(text = paste0("!is.na(",i,")")))), fill = i) +
      geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
      ylab("Abundance (read count)") + scale_fill_manual(values = temp_col_pallete, na.value="lightgrey") + 
      scale_color_manual(values = temp_col_pallete, na.value="lightgrey") +
      theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
  } else {
    bar_plot <- plot_bar(phylo4, fill = i) +
      geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
      ylab("Abundance (read count)") + scale_fill_manual(values = temp_col_pallete, na.value="lightgrey") + 
      scale_color_manual(values = temp_col_pallete, na.value="lightgrey") +
      theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
  }
  bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = sample_order$V1)
  bar_plot_legend <- get_legend(bar_plot)
  bar_plot <- bar_plot + theme(legend.position='none')
  pdf(file=paste0('barplot_readcount_filtsamples_filtLT',filter_percent,'PERCtaxa_',i,'_legend.pdf', sep = ""), width = 11, height = 6)
  print(as_ggplot(bar_plot_legend))
  dev.off()
  pdf(file=paste0('barplot_readcount_filtsamples_filtLT',filter_percent,'PERCtaxa_',i,'.pdf', sep = ""), width = 11, height = 8.5)
  print(bar_plot)
  dev.off()
}
}, error=function(e){})

###Relative abundance
setwd(paste0(as.character(args[1]), "/02_Barcharts/relative_abundance", sep = ""))

tryCatch({
for (i in c("Phylum")) {
  if (NAremoveFlag == TRUE) {
    standf = function(x) 100*(x / sum(x))
    bar_plot <- plot_bar(transform_sample_counts(subset_taxa(phylo3, eval(parse(text = paste0("!is.na(",i,")")))), standf), fill = i) +
      geom_bar(aes(fill=.data[[i]]), stat="identity", position="stack") +
      ylab("Relative Abundance (%)") +
      theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
  } else {
    bar_plot <- plot_bar(phylo3, fill = i) +
      geom_bar(aes(fill=.data[[i]]), stat="identity", position="stack") +
      ylab("Relative Abundance (%)") +
      theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
  }
  bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = sample_order$V1)
  bar_plot <- bar_plot + theme(legend.position='none')
  pdf(file=paste0('barplot_relabund_filtsamples_alltaxa_',i,'_ASVsShown.pdf', sep = ""), width = 11, height = 8.5)
  print(bar_plot)
  dev.off()
}
}, error=function(e){})

tryCatch({
for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
  if (NAremoveFlag == TRUE) {
    standf = function(x) 100*(x / sum(x))
    bar_plot <- plot_bar(transform_sample_counts(subset_taxa(phylo3, eval(parse(text = paste0("!is.na(",i,")")))), standf), fill = i) +
      geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
      ylab("Relative Abundance (%)") +
      theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
  } else {
    bar_plot <- plot_bar(phylo3, fill = i) +
      geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
      ylab("Relative Abundance (%)") +
      theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
  }
  bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = sample_order$V1)
  bar_plot_legend <- get_legend(bar_plot)
  bar_plot <- bar_plot + theme(legend.position='none')
  legend_call <- eval(parse(text = paste0("base_",i,"_count")))
  w_choice <- 6.4 + ((legend_call*2.3)/20)
  pdf(file=paste0('barplot_relabund_filtsamples_alltaxa_',i,'_legend.pdf', sep = ""), width = w_choice, height = 6)
  print(as_ggplot(bar_plot_legend))
  dev.off()
  pdf(file=paste0('barplot_relabund_filtsamples_alltaxa_',i,'.pdf', sep = ""), width = 11, height = 8.5)
  print(bar_plot)
  dev.off()
}
}, error=function(e){})

tryCatch({
for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
  count_call <- eval(parse(text = paste0(i,"_count")))
  NA_call <- eval(parse(text = paste0(i,"_NA")))
  temp_col_pallete <- c(gg_color(count_call-NA_call-1), "grey")
  if (NAremoveFlag == TRUE) {
    standf = function(x) 100*(x / sum(x))
    bar_plot <- plot_bar(transform_sample_counts(subset_taxa(phylo5, eval(parse(text = paste0("!is.na(",i,")")))), standf), fill = i) +
      geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
      ylab("Relative Abundance (%)") + scale_fill_manual(values = temp_col_pallete, na.value="lightgrey") + 
      scale_color_manual(values = temp_col_pallete, na.value="lightgrey") +
      theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
  } else {
    bar_plot <- plot_bar(phylo5, fill = i) +
      geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
      ylab("Relative Abundance (%)") + scale_fill_manual(values = temp_col_pallete, na.value="lightgrey") + 
      scale_color_manual(values = temp_col_pallete, na.value="lightgrey") +
      theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
  }
  bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = sample_order$V1)
  bar_plot_legend <- get_legend(bar_plot)
  bar_plot <- bar_plot + theme(legend.position='none')
  pdf(file=paste0('barplot_relabund_filtsamples_filtLT',filter_percent,'PERCtaxa_',i,'_legend.pdf', sep = ""), width = 11, height = 6)
  print(as_ggplot(bar_plot_legend))
  dev.off()
  pdf(file=paste0('barplot_relabund_filtsamples_filtLT',filter_percent,'PERCtaxa_',i,'.pdf', sep = ""), width = 11, height = 8.5)
  print(bar_plot)
  dev.off()
}
}, error=function(e){})

#GROUPED
tryCatch({
if (replicateFlag == TRUE) {
  for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    if (NAremoveFlag == TRUE) {
      standf = function(x) 100*(x / sum(x))
      bar_plot <- plot_bar(transform_sample_counts(subset_taxa(phylo7a, eval(parse(text = paste0("!is.na(",i,")")))), standf), fill = i) +
        geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
        ylab("Relative Abundance (%)") +
        theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
    } else {
      bar_plot <- plot_bar(phylo7a, fill = i) +
        geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
        ylab("Relative Abundance (%)") +
        theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
    }
    bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = replicateorder$`sample_metadata$replicates`)
    bar_plot_legend <- get_legend(bar_plot)
    bar_plot <- bar_plot + theme(legend.position='none')
    legend_call <- eval(parse(text = paste0("base_",i,"_count")))
    w_choice <- 6.4 + ((legend_call*2.3)/20)
    pdf(file=paste0('barplot_relabund_replicateGroupedSamples_alltaxa_',i,'_legend.pdf', sep = ""), width = w_choice, height = 6)
    print(as_ggplot(bar_plot_legend))
    dev.off()
    pdf(file=paste0('barplot_relabund_replicateGroupedSamples_alltaxa_',i,'.pdf', sep = ""), width = 11, height = 8.5)
    print(bar_plot)
    dev.off()
  }
}
}, error=function(e){})

tryCatch({
if (sitelabelFlag == TRUE) {
  for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    if (NAremoveFlag == TRUE) {
      standf = function(x) 100*(x / sum(x))
      bar_plot <- plot_bar(transform_sample_counts(subset_taxa(phylo7b, eval(parse(text = paste0("!is.na(",i,")")))), standf), fill = i) +
        geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
        ylab("Relative Abundance (%)") +
        theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
    } else {
      bar_plot <- plot_bar(phylo7b, fill = i) +
        geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
        ylab("Relative Abundance (%)") +
        theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
    }
    bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = siteorder$`sample_metadata$sites`)
    bar_plot_legend <- get_legend(bar_plot)
    bar_plot <- bar_plot + theme(legend.position='none')
    legend_call <- eval(parse(text = paste0("base_",i,"_count")))
    w_choice <- 6.4 + ((legend_call*2.3)/20)
    pdf(file=paste0('barplot_relabund_siteGroupedSamples_alltaxa_',i,'_legend.pdf', sep = ""), width = w_choice, height = 6)
    print(as_ggplot(bar_plot_legend))
    dev.off()
    pdf(file=paste0('barplot_relabund_siteGroupedSamples_alltaxa_',i,'.pdf', sep = ""), width = 11, height = 8.5)
    print(bar_plot)
    dev.off()
  }
}
}, error=function(e){})

tryCatch({
if (replicateFlag == TRUE) {
  for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    count_call <- eval(parse(text = paste0(i,"_count")))
    NA_call <- eval(parse(text = paste0(i,"_NA")))
    temp_col_pallete <- c(gg_color(count_call-NA_call-1), "grey")
    if (NAremoveFlag == TRUE) {
      standf = function(x) 100*(x / sum(x))
      bar_plot <- plot_bar(transform_sample_counts(subset_taxa(phylo9a, eval(parse(text = paste0("!is.na(",i,")")))), standf), fill = i) +
        geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
        ylab("Relative Abundance (%)") + scale_fill_manual(values = temp_col_pallete, na.value="lightgrey") + 
        scale_color_manual(values = temp_col_pallete, na.value="lightgrey") +
        theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
    } else {
      bar_plot <- plot_bar(phylo9a, fill = i) +
        geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
        ylab("Relative Abundance (%)") + scale_fill_manual(values = temp_col_pallete, na.value="lightgrey") + 
        scale_color_manual(values = temp_col_pallete, na.value="lightgrey") +
        theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
    }
    bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = replicateorder$`sample_metadata$replicates`)
    bar_plot_legend <- get_legend(bar_plot)
    bar_plot <- bar_plot + theme(legend.position='none')
    pdf(file=paste0('barplot_relabund_replicateGroupedSamples_filtLT',filter_percent,'PERCtaxa_',i,'_legend.pdf', sep = ""), width = 11, height = 6)
    print(as_ggplot(bar_plot_legend))
    dev.off()
    pdf(file=paste0('barplot_relabund_replicateGroupedSamples_filtLT',filter_percent,'PERCtaxa_',i,'.pdf', sep = ""), width = 11, height = 8.5)
    print(bar_plot)
    dev.off()
  }
}
}, error=function(e){})

tryCatch({
if (sitelabelFlag == TRUE) {
  for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    count_call <- eval(parse(text = paste0(i,"_count")))
    NA_call <- eval(parse(text = paste0(i,"_NA")))
    temp_col_pallete <- c(gg_color(count_call-NA_call-1), "grey")
    if (NAremoveFlag == TRUE) {
      standf = function(x) 100*(x / sum(x))
      bar_plot <- plot_bar(transform_sample_counts(subset_taxa(phylo9b, eval(parse(text = paste0("!is.na(",i,")")))), standf), fill = i) +
        geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
        ylab("Relative Abundance (%)") + scale_fill_manual(values = temp_col_pallete, na.value="lightgrey") + 
        scale_color_manual(values = temp_col_pallete, na.value="lightgrey") +
        theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
    } else {
      bar_plot <- plot_bar(phylo9b, fill = i) +
        geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
        ylab("Relative Abundance (%)") + scale_fill_manual(values = temp_col_pallete, na.value="lightgrey") + 
        scale_color_manual(values = temp_col_pallete, na.value="lightgrey") +
        theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
    }
    bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = siteorder$`sample_metadata$sites`)
    bar_plot_legend <- get_legend(bar_plot)
    bar_plot <- bar_plot + theme(legend.position='none')
    pdf(file=paste0('barplot_relabund_siteGroupedSamples_filtLT',filter_percent,'PERCtaxa_',i,'_legend.pdf', sep = ""), width = 11, height = 6)
    print(as_ggplot(bar_plot_legend))
    dev.off()
    pdf(file=paste0('barplot_relabund_siteGroupedSamples_filtLT',filter_percent,'PERCtaxa_',i,'.pdf', sep = ""), width = 11, height = 8.5)
    print(bar_plot)
    dev.off()
  }
}
}, error=function(e){})
}, error=function(e){})
########################
#
# HEATMAPS
#
########################
setwd(paste0(as.character(args[1]), "/03_Heatmaps/ASV_based", sep = ""))

tryCatch({
heatmap <- plot_heatmap(phylo3, method = "NMDS", distance = "jaccard", 
             low = "ghostwhite", 
             high = "red", 
             na.value = "white") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
tryCatch({
heatmap$scales$scales[[1]]$name <- "ASV"
heatmap$scales$scales[[2]]$name <- "Relative Abundance (%)"
}, error=function(e){heatmap$scales$scales[[1]]$name <- "Relative Abundance (%)" })
pdf(file='heatmap_ASVs_relabund_filtsamples_alltaxa_clustSamples.pdf', width = 11, height = 8.5)
print(heatmap)
dev.off()
heatmap$data$Sample <- factor(heatmap$data$Sample, levels = sample_order$V1)
pdf(file='heatmap_ASVs_relabund_filtsamples_alltaxa_orderedSamples.pdf', width = 11, height = 8.5)
print(heatmap)
dev.off()
}, error=function(e){})

if (replicateFlag == TRUE) {
  tryCatch({
  heatmap <- plot_heatmap(phylo7a, method = "NMDS", distance = "jaccard", 
                          low = "ghostwhite", 
                          high = "red", 
                          na.value = "white") +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
  tryCatch({
  heatmap$scales$scales[[1]]$name <- "ASV"
  heatmap$scales$scales[[2]]$name <- "Relative Abundance (%)"
  }, error=function(e){heatmap$scales$scales[[1]]$name <- "Relative Abundance (%)" })
  pdf(file='heatmap_ASVs_relabund_replicateGroupedSamples_alltaxa_clustReplicates.pdf', width = 11, height = 8.5)
  print(heatmap)
  dev.off()
  heatmap$data$Sample <- factor(heatmap$data$Sample, levels = replicateorder$`sample_metadata$replicates`)
  pdf(file='heatmap_ASVs_relabund_replicateGroupedSamples_alltaxa_orderedReplicates.pdf', width = 11, height = 8.5)
  print(heatmap)
  dev.off()
  }, error=function(e){})
}


if (sitelabelFlag == TRUE) {
  tryCatch({
  heatmap <- plot_heatmap(phylo7b, method = "NMDS", distance = "jaccard", 
                          low = "ghostwhite", 
                          high = "red", 
                          na.value = "white") +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
  tryCatch({
  heatmap$scales$scales[[1]]$name <- "ASV"
  heatmap$scales$scales[[2]]$name <- "Relative Abundance (%)"
  }, error=function(e){heatmap$scales$scales[[1]]$name <- "Relative Abundance (%)" })
  pdf(file='heatmap_ASVs_relabund_siteGroupedSamples_alltaxa_clustSites.pdf', width = 11, height = 8.5)
  print(heatmap)
  dev.off()
  heatmap$data$Sample <- factor(heatmap$data$Sample, levels = siteorder$`sample_metadata$sites`)
  pdf(file='heatmap_ASVs_relabund_siteGroupedSamples_alltaxa_orderedSites.pdf', width = 11, height = 8.5)
  print(heatmap)
  dev.off()
  }, error=function(e){})
}

setwd(paste0(as.character(args[1]), "/03_Heatmaps/Taxonomy_merge_based", sep = ""))

tryCatch({
heatmap <- plot_heatmap(phylo11, method = "NMDS", distance = "jaccard", 
                        low = "ghostwhite", 
                        high = "red", 
                        na.value = "white") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
tryCatch({
heatmap$scales$scales[[1]]$name <- "ASV"
heatmap$scales$scales[[2]]$name <- "Relative Abundance (%)"
}, error=function(e){heatmap$scales$scales[[1]]$name <- "Relative Abundance (%)" })
pdf(file='heatmap_Taxa_relabund_filtsamples_alltaxa_clustSamples.pdf', width = 11, height = 8.5)
print(heatmap)
dev.off()
heatmap$data$Sample <- factor(heatmap$data$Sample, levels = sample_order$V1)
pdf(file='heatmap_Taxa_relabund_filtsamples_alltaxa_orderedSamples.pdf', width = 11, height = 8.5)
print(heatmap)
dev.off()
}, error=function(e){})

if (replicateFlag == TRUE) {
  tryCatch({
  heatmap <- plot_heatmap(phylo15a, method = "NMDS", distance = "jaccard", 
                          low = "ghostwhite", 
                          high = "red", 
                          na.value = "white") +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
  tryCatch({
  heatmap$scales$scales[[1]]$name <- "ASV"
  heatmap$scales$scales[[2]]$name <- "Relative Abundance (%)"
  }, error=function(e){heatmap$scales$scales[[1]]$name <- "Relative Abundance (%)" })
  pdf(file='heatmap_Taxa_relabund_replicateGroupedSamples_alltaxa_clustReplicates.pdf', width = 11, height = 8.5)
  print(heatmap)
  dev.off()
  heatmap$data$Sample <- factor(heatmap$data$Sample, levels = replicateorder$`sample_metadata$replicates`)
  pdf(file='heatmap_Taxa_relabund_replicateGroupedSamples_alltaxa_orderedReplicates.pdf', width = 11, height = 8.5)
  print(heatmap)
  dev.off()
  }, error=function(e){})
}


if (sitelabelFlag == TRUE) {
  tryCatch({
  heatmap <- plot_heatmap(phylo15b, method = "NMDS", distance = "jaccard", 
                          low = "ghostwhite", 
                          high = "red", 
                          na.value = "white") +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
  tryCatch({
  heatmap$scales$scales[[1]]$name <- "ASV"
  heatmap$scales$scales[[2]]$name <- "Relative Abundance (%)"
  }, error=function(e){heatmap$scales$scales[[1]]$name <- "Relative Abundance (%)" })
  pdf(file='heatmap_Taxa_relabund_siteGroupedSamples_alltaxa_clustSites.pdf', width = 11, height = 8.5)
  print(heatmap)
  dev.off()
  heatmap$data$Sample <- factor(heatmap$data$Sample, levels = siteorder$`sample_metadata$sites`)
  pdf(file='heatmap_Taxa_relabund_siteGroupedSamples_alltaxa_orderedSites.pdf', width = 11, height = 8.5)
  print(heatmap)
  dev.off()
  }, error=function(e){})
}

########################
#
# ALPHA DIVERSITY
# (from normalized dataset only)
#
########################
setwd(paste0(as.character(args[1]), "/04_Alpha_Diversity/ASV_based", sep = ""))

tryCatch({
total = median(sample_sums(phylo2))
standf = function(x, t=total) round(t * (x / sum(x)))

tryCatch({
alpha <- plot_richness(transform_sample_counts(phylo2, standf))
out_alpha <- as.data.frame(alpha$data)
write.table(out_alpha, "alpha_diversity_ASVbased_R.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}, error=function(e){})

tryCatch({
alpha_human <- estimate_richness(transform_sample_counts(phylo2, standf))
for (i in c("Chao1", "Shannon")) {
  if (replicateFlag == TRUE) {
    message(paste0("Wilcox Rank Sum Test: Significant difference between Replicates based on ",i," for phylo2 (ASV-based)"))
    print(pairwise.wilcox.test(alpha_human[[i]], sample_data(phylo2)$replicates))
  }
  if (sitelabelFlag == TRUE) {
    message(paste0("Wilcox Rank Sum Test: Significant difference between Sites based on ",i," for phylo2 (ASV-based)"))
    print(pairwise.wilcox.test(alpha_human[[i]], sample_data(phylo2)$sites))
  }
  if (groupingFlag == TRUE) {
    for (j in 1:numberofGroups) {
      message(paste0("Wilcox Rank Sum Test: Significant difference between Group ",j," based on ",i," for phylo2 (ASV-based)"))
      print(pairwise.wilcox.test(alpha_human[[i]], sample_data(phylo2)[[paste0("group",j)]]))
    }
  }
}
alpha_human <- tibble::rownames_to_column(alpha_human, "Sample")
write.table(alpha_human, "alpha_diversity_ASVbased_human.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}, error=function(e){})

tryCatch({
if (replicateFlag == TRUE) {
  total = median(sample_sums(phylo2))
  alpha <- plot_richness(transform_sample_counts(phylo2, standf), nrow = NULL, x = "replicates") +
    xlab("Replicates")
  alpha$data$replicates <- factor(alpha$data$replicates, levels = replicateorder$`sample_metadata$replicates`)
  pdf(file='alphaDiversity_ASVs_normalized_replicateClusteredSamples_alltaxa_allMeasures.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
  
  alpha <- plot_richness(transform_sample_counts(phylo2, standf), measures = c("Chao1","Shannon"), x = "replicates") +
    xlab("Replicates")
  alpha$data$replicates <- factor(alpha$data$replicates, levels = replicateorder$`sample_metadata$replicates`)
  pdf(file='alphaDiversity_ASVs_normalized_replicateClusteredSamples_alltaxa_Chao1Shannon.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
  
  total = median(sample_sums(phylo7a))
  alpha <- plot_richness(transform_sample_counts(phylo7a, standf), nrow = NULL) +
    xlab("Replicates")
  alpha$data$samples <- factor(alpha$data$samples, levels = replicateorder$`sample_metadata$replicates`)
  pdf(file='alphaDiversity_ASVs_normalized_replicateGroupedSamples_alltaxa_allMeasures.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
  
  alpha <- plot_richness(transform_sample_counts(phylo7a, standf), measures = c("Chao1","Shannon")) +
    xlab("Replicates")
  alpha$data$samples <- factor(alpha$data$samples, levels = replicateorder$`sample_metadata$replicates`)
  pdf(file='alphaDiversity_ASVs_normalized_replicateGroupedSamples_alltaxa_Chao1Shannon.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
}
}, error=function(e){})

tryCatch({
if (sitelabelFlag == TRUE) {
  total = median(sample_sums(phylo2))
  alpha <- plot_richness(transform_sample_counts(phylo2, standf), nrow = NULL, x = "sites") +
    xlab("Sites")
  alpha$data$sites <- factor(alpha$data$sites, levels = siteorder$`sample_metadata$sites`)
  pdf(file='alphaDiversity_ASVs_normalized_siteClusteredSamples_alltaxa_allMeasures.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
  
  alpha <- plot_richness(transform_sample_counts(phylo2, standf), measures = c("Chao1","Shannon"), x = "sites") +
    xlab("Sites")
  alpha$data$sites <- factor(alpha$data$sites, levels = siteorder$`sample_metadata$sites`)
  pdf(file='alphaDiversity_ASVs_normalized_siteClusteredSamples_alltaxa_Chao1Shannon.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
  
  total = median(sample_sums(phylo7b))
  alpha <- plot_richness(transform_sample_counts(phylo7b, standf), nrow = NULL) +
    xlab("Sites")
  alpha$data$samples <- factor(alpha$data$samples, levels = siteorder$`sample_metadata$sites`)
  pdf(file='alphaDiversity_ASVs_normalized_siteGroupedSamples_alltaxa_allMeasures.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
  
  alpha <- plot_richness(transform_sample_counts(phylo7b, standf), measures = c("Chao1","Shannon")) +
    xlab("Sites")
  alpha$data$samples <- factor(alpha$data$samples, levels = siteorder$`sample_metadata$sites`)
  pdf(file='alphaDiversity_ASVs_normalized_siteGroupedSamples_alltaxa_Chao1Shannon.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
}
}, error=function(e){})

tryCatch({
total = median(sample_sums(phylo2))
alpha <- plot_richness(transform_sample_counts(phylo2, standf), nrow = NULL) +
  xlab("Samples")
alpha$data$samples <- factor(alpha$data$samples, levels = sample_order$V1)
pdf(file='alphaDiversity_ASVs_normalized_filtsamples_alltaxa_allMeasures.pdf', width = 15, height = 8.5)
print(alpha)
dev.off()

alpha <- plot_richness(transform_sample_counts(phylo2, standf), measures = c("Chao1","Shannon")) +
  xlab("Samples")
alpha$data$samples <- factor(alpha$data$samples, levels = sample_order$V1)
pdf(file='alphaDiversity_ASVs_normalized_filtsamples_alltaxa_Chao1Shannon.pdf', width = 15, height = 8.5)
print(alpha)
dev.off()
}, error=function(e){})


if (groupingFlag == TRUE) {
  for (i in 1:numberofGroups) {
    tryCatch({
    if (replicateFlag == TRUE) {
      total = median(sample_sums(phylo2))
      alpha <- plot_richness(transform_sample_counts(phylo2, standf), nrow = NULL, x = "replicates", color = paste0("group",i)) +
        xlab("Replicates")
      alpha$data$replicates <- factor(alpha$data$replicates, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_ASVs_normalized_replicateClusteredSamples_alltaxa_allMeasures_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(phylo2, standf), measures = c("Chao1","Shannon"), x = "replicates", color = paste0("group",i)) +
        xlab("Replicates")
      alpha$data$replicates <- factor(alpha$data$replicates, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_ASVs_normalized_replicateClusteredSamples_alltaxa_Chao1Shannon_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      total = median(sample_sums(phylo7a))
      alpha <- plot_richness(transform_sample_counts(phylo7a, standf), nrow = NULL, color = paste0("group",i)) +
        xlab("Replicates")
      alpha$data$samples <- factor(alpha$data$samples, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_ASVs_normalized_replicateGroupedSamples_alltaxa_allMeasures_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(phylo7a, standf), measures = c("Chao1","Shannon"), color = paste0("group",i)) +
        xlab("Replicates")
      alpha$data$samples <- factor(alpha$data$samples, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_ASVs_normalized_replicateGroupedSamples_alltaxa_Chao1Shannon_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    if (sitelabelFlag == TRUE) {
      total = median(sample_sums(phylo2))
      alpha <- plot_richness(transform_sample_counts(phylo2, standf), nrow = NULL, x = "sites", color = paste0("group",i)) +
        xlab("Sites")
      alpha$data$sites <- factor(alpha$data$sites, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_ASVs_normalized_siteClusteredSamples_alltaxa_allMeasures_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(phylo2, standf), measures = c("Chao1","Shannon"), x = "sites", color = paste0("group",i)) +
        xlab("Sites")
      alpha$data$sites <- factor(alpha$data$sites, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_ASVs_normalized_siteClusteredSamples_alltaxa_Chao1Shannon_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      total = median(sample_sums(phylo7b))
      alpha <- plot_richness(transform_sample_counts(phylo7b, standf), nrow = NULL, color = paste0("group",i)) +
        xlab("Sites")
      alpha$data$samples <- factor(alpha$data$samples, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_ASVs_normalized_siteGroupedSamples_alltaxa_allMeasures_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(phylo7b, standf), measures = c("Chao1","Shannon"), color = paste0("group",i)) +
        xlab("Sites")
      alpha$data$samples <- factor(alpha$data$samples, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_ASVs_normalized_siteGroupedSamples_alltaxa_Chao1Shannon_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    total = median(sample_sums(phylo2))
    alpha <- plot_richness(transform_sample_counts(phylo2, standf), nrow = NULL, color = paste0("group",i), sortby = "Chao1") +
      xlab("Samples")
    pdf(file=paste0('alphaDiversity_ASVs_normalized_filtsamples_alltaxa_allMeasures_group',i,'Colored_sortonChao1.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha <- plot_richness(transform_sample_counts(phylo2, standf), nrow = NULL, color = paste0("group",i), sortby = "Shannon") +
      xlab("Samples")
    pdf(file=paste0('alphaDiversity_ASVs_normalized_filtsamples_alltaxa_allMeasures_group',i,'Colored_sortonShannon.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha$data$samples <- factor(alpha$data$samples, levels = sample_order$V1)
    pdf(file=paste0('alphaDiversity_ASVs_normalized_filtsamples_alltaxa_allMeasures_group',i,'Colored.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    
    alpha <- plot_richness(transform_sample_counts(phylo2, standf), nrow = NULL, color = paste0("group",i), measures = c("Chao1","Shannon"), sortby = "Chao1") +
      xlab("Samples")
    pdf(file=paste0('alphaDiversity_ASVs_normalized_filtsamples_alltaxa_Chao1Shannon_group',i,'Colored_sortonChao1.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha <- plot_richness(transform_sample_counts(phylo2, standf), nrow = NULL, color = paste0("group",i), measures = c("Chao1","Shannon"), sortby = "Shannon") +
      xlab("Samples")
    pdf(file=paste0('alphaDiversity_ASVs_normalized_filtsamples_alltaxa_Chao1Shannon_group',i,'Colored_sortonShannon.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha$data$samples <- factor(alpha$data$samples, levels = sample_order$V1)
    pdf(file=paste0('alphaDiversity_ASVs_normalized_filtsamples_alltaxa_Chao1Shannon_group',i,'Colored.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    }, error=function(e){})
  }
}

if (chemDataFlag == TRUE) {
  for (i in chem_headers$V1) {
    tryCatch({
    if (replicateFlag == TRUE) {
      classgroup <- class(sample_data(phylo2)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      total = median(sample_sums(phylo2))
      alpha <- plot_richness(transform_sample_counts(phylo2, standf), nrow = NULL, x = "replicates", color = paste0(i)) +
        xlab("Replicates") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$replicates <- factor(alpha$data$replicates, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_ASVs_normalized_replicateClusteredSamples_alltaxa_allMeasures_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(phylo2, standf), measures = c("Chao1","Shannon"), x = "replicates", color = paste0(i)) +
        xlab("Replicates") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$replicates <- factor(alpha$data$replicates, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_ASVs_normalized_replicateClusteredSamples_alltaxa_Chao1Shannon_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      classgroup <- class(sample_data(phylo7a)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      total = median(sample_sums(phylo7a))
      alpha <- plot_richness(transform_sample_counts(phylo7a, standf), nrow = NULL, color = paste0(i)) +
        xlab("Replicates") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$samples <- factor(alpha$data$samples, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_ASVs_normalized_replicateGroupedSamples_alltaxa_allMeasures_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(phylo7a, standf), measures = c("Chao1","Shannon"), color = paste0(i)) +
        xlab("Replicates") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$samples <- factor(alpha$data$samples, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_ASVs_normalized_replicateGroupedSamples_alltaxa_Chao1Shannon_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    if (sitelabelFlag == TRUE) {
      classgroup <- class(sample_data(phylo2)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      total = median(sample_sums(phylo2))
      alpha <- plot_richness(transform_sample_counts(phylo2, standf), nrow = NULL, x = "sites", color = paste0(i)) +
        xlab("Sites") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$sites <- factor(alpha$data$sites, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_ASVs_normalized_siteClusteredSamples_alltaxa_allMeasures_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(phylo2, standf), measures = c("Chao1","Shannon"), x = "sites", color = paste0(i)) +
        xlab("Sites") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$sites <- factor(alpha$data$sites, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_ASVs_normalized_siteClusteredSamples_alltaxa_Chao1Shannon_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      classgroup <- class(sample_data(phylo7b)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      total = median(sample_sums(phylo7b))
      alpha <- plot_richness(transform_sample_counts(phylo7b, standf), nrow = NULL, color = paste0(i)) +
        xlab("Sites") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$samples <- factor(alpha$data$samples, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_ASVs_normalized_siteGroupedSamples_alltaxa_allMeasures_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(phylo7b, standf), measures = c("Chao1","Shannon"), color = paste0(i)) +
        xlab("Sites") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$samples <- factor(alpha$data$samples, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_ASVs_normalized_siteGroupedSamples_alltaxa_Chao1Shannon_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    classgroup <- class(sample_data(phylo2)[[i]])
    if (classgroup != "integer" && classgroup != "numeric") {
      viridis_scale <- 'scale_color_viridis_d(option="C")'
    } else {
      viridis_scale <- 'scale_color_viridis_c(option="C")'
    }
    
    total = median(sample_sums(phylo2))
    alpha <- plot_richness(transform_sample_counts(phylo2, standf), nrow = NULL, color = paste0(i), sortby = "Chao1") +
      xlab("Samples") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
    pdf(file=paste0('alphaDiversity_ASVs_normalized_filtsamples_alltaxa_allMeasures_chem_',i,'_Colored_sortonChao1.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha <- plot_richness(transform_sample_counts(phylo2, standf), nrow = NULL, color = paste0(i), sortby = "Shannon") +
      xlab("Samples") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
    pdf(file=paste0('alphaDiversity_ASVs_normalized_filtsamples_alltaxa_allMeasures_chem_',i,'_Colored_sortonShannon.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha$data$samples <- factor(alpha$data$samples, levels = sample_order$V1)
    pdf(file=paste0('alphaDiversity_ASVs_normalized_filtsamples_alltaxa_allMeasures_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    
    alpha <- plot_richness(transform_sample_counts(phylo2, standf), nrow = NULL, color = paste0(i), measures = c("Chao1","Shannon"), sortby = "Chao1") +
      xlab("Samples") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
    pdf(file=paste0('alphaDiversity_ASVs_normalized_filtsamples_alltaxa_Chao1Shannon_chem_',i,'_Colored_sortonChao1.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha <- plot_richness(transform_sample_counts(phylo2, standf), nrow = NULL, color = paste0(i), measures = c("Chao1","Shannon"), sortby = "Shannon") +
      xlab("Samples") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
    pdf(file=paste0('alphaDiversity_ASVs_normalized_filtsamples_alltaxa_Chao1Shannon_chem_',i,'_Colored_sortonShannon.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha$data$samples <- factor(alpha$data$samples, levels = sample_order$V1)
    pdf(file=paste0('alphaDiversity_ASVs_normalized_filtsamples_alltaxa_Chao1Shannon_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    }, error=function(e){})
  }
}

###############
# Taxa-based
###############

setwd(paste0(as.character(args[1]), "/04_Alpha_Diversity/Taxonomy_merge_based", sep = ""))

total = median(sample_sums(phylo10))
standf = function(x, t=total) round(t * (x / sum(x)))

tryCatch({
alpha <- plot_richness(transform_sample_counts(phylo10, standf))
out_alpha <- as.data.frame(alpha$data)
write.table(out_alpha, "alpha_diversity_TAXAbased_R.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}, error=function(e){})

tryCatch({
alpha_human <- estimate_richness(transform_sample_counts(phylo10, standf))
for (i in c("Chao1", "Shannon")) {
  if (replicateFlag == TRUE) {
    message(paste0("Wilcox Rank Sum Test: Significant difference between Replicates based on ",i," for phylo10 (TAXA-based)"))
    print(pairwise.wilcox.test(alpha_human[[i]], sample_data(phylo10)$replicates))
  }
  if (sitelabelFlag == TRUE) {
    message(paste0("Wilcox Rank Sum Test: Significant difference between Sites based on ",i," for phylo10 (TAXA-based)"))
    print(pairwise.wilcox.test(alpha_human[[i]], sample_data(phylo10)$sites))
  }
  if (groupingFlag == TRUE) {
    for (j in 1:numberofGroups) {
      message(paste0("Wilcox Rank Sum Test: Significant difference between Group ",j," based on ",i," for phylo10 (TAXA-based)"))
      print(pairwise.wilcox.test(alpha_human[[i]], sample_data(phylo10)[[paste0("group",j)]]))
    }
  }
}
alpha_human <- tibble::rownames_to_column(alpha_human, "Sample")
write.table(alpha_human, "alpha_diversity_TAXAbased_human.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}, error=function(e){})

tryCatch({
alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf))
out_alpha <- as.data.frame(alpha$data)
write.table(out_alpha, "alpha_diversity_TAXAbased_ONLYINCLUDINGTAXACOMPLETETO_Species_R.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
alpha_human <- estimate_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf))
alpha_human <- tibble::rownames_to_column(alpha_human, "Sample")
write.table(alpha_human, "alpha_diversity_TAXAbased_ONLYINCLUDINGTAXACOMPLETETO_Species_human.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}, error=function(e){})

tryCatch({
alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Genus)), standf))
out_alpha <- as.data.frame(alpha$data)
write.table(out_alpha, "alpha_diversity_TAXAbased_ONLYINCLUDINGTAXACOMPLETETO_Genus_R.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
alpha_human <- estimate_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Genus)), standf))
alpha_human <- tibble::rownames_to_column(alpha_human, "Sample")
write.table(alpha_human, "alpha_diversity_TAXAbased_ONLYINCLUDINGTAXACOMPLETETO_Genus_human.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}, error=function(e){})

tryCatch({
alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Family)), standf))
out_alpha <- as.data.frame(alpha$data)
write.table(out_alpha, "alpha_diversity_TAXAbased_ONLYINCLUDINGTAXACOMPLETETO_Family_R.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
alpha_human <- estimate_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Family)), standf))
alpha_human <- tibble::rownames_to_column(alpha_human, "Sample")
write.table(alpha_human, "alpha_diversity_TAXAbased_ONLYINCLUDINGTAXACOMPLETETO_Family_human.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}, error=function(e){})

tryCatch({
alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Order)), standf))
out_alpha <- as.data.frame(alpha$data)
write.table(out_alpha, "alpha_diversity_TAXAbased_ONLYINCLUDINGTAXACOMPLETETO_Order_R.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
alpha_human <- estimate_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Order)), standf))
alpha_human <- tibble::rownames_to_column(alpha_human, "Sample")
write.table(alpha_human, "alpha_diversity_TAXAbased_ONLYINCLUDINGTAXACOMPLETETO_Order_human.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}, error=function(e){})

tryCatch({
alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Class)), standf))
out_alpha <- as.data.frame(alpha$data)
write.table(out_alpha, "alpha_diversity_TAXAbased_ONLYINCLUDINGTAXACOMPLETETO_Class_R.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
alpha_human <- estimate_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Class)), standf))
alpha_human <- tibble::rownames_to_column(alpha_human, "Sample")
write.table(alpha_human, "alpha_diversity_TAXAbased_ONLYINCLUDINGTAXACOMPLETETO_Class_human.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}, error=function(e){})

tryCatch({
alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Phylum)), standf))
out_alpha <- as.data.frame(alpha$data)
write.table(out_alpha, "alpha_diversity_TAXAbased_ONLYINCLUDINGTAXACOMPLETETO_Phylum_R.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
alpha_human <- estimate_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Phylum)), standf))
alpha_human <- tibble::rownames_to_column(alpha_human, "Sample")
write.table(alpha_human, "alpha_diversity_TAXAbased_ONLYINCLUDINGTAXACOMPLETETO_Phylum_human.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}, error=function(e){})

tryCatch({
if (replicateFlag == TRUE) {
  total = median(sample_sums(phylo10))
  alpha <- plot_richness(transform_sample_counts(phylo10, standf), nrow = NULL, x = "replicates") +
    xlab("Replicates")
  alpha$data$replicates <- factor(alpha$data$replicates, levels = replicateorder$`sample_metadata$replicates`)
  pdf(file='alphaDiversity_Taxa_normalized_replicateClusteredSamples_alltaxa_allMeasures.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
  
  alpha <- plot_richness(transform_sample_counts(phylo10, standf), measures = c("Chao1","Shannon"), x = "replicates") +
    xlab("Replicates")
  alpha$data$replicates <- factor(alpha$data$replicates, levels = replicateorder$`sample_metadata$replicates`)
  pdf(file='alphaDiversity_Taxa_normalized_replicateClusteredSamples_alltaxa_Chao1Shannon.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
  
  total = median(sample_sums(phylo15a))
  alpha <- plot_richness(transform_sample_counts(phylo15a, standf), nrow = NULL) +
    xlab("Replicates")
  alpha$data$samples <- factor(alpha$data$samples, levels = replicateorder$`sample_metadata$replicates`)
  pdf(file='alphaDiversity_Taxa_normalized_replicateGroupedSamples_alltaxa_allMeasures.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
  
  alpha <- plot_richness(transform_sample_counts(phylo15a, standf), measures = c("Chao1","Shannon")) +
    xlab("Replicates")
  alpha$data$samples <- factor(alpha$data$samples, levels = replicateorder$`sample_metadata$replicates`)
  pdf(file='alphaDiversity_Taxa_normalized_replicateGroupedSamples_alltaxa_Chao1Shannon.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
}
}, error=function(e){})

tryCatch({
if (sitelabelFlag == TRUE) {
  total = median(sample_sums(phylo10))
  alpha <- plot_richness(transform_sample_counts(phylo10, standf), nrow = NULL, x = "sites") +
    xlab("Sites")
  alpha$data$sites <- factor(alpha$data$sites, levels = siteorder$`sample_metadata$sites`)
  pdf(file='alphaDiversity_Taxa_normalized_siteClusteredSamples_alltaxa_allMeasures.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
  
  alpha <- plot_richness(transform_sample_counts(phylo10, standf), measures = c("Chao1","Shannon"), x = "sites") +
    xlab("Sites")
  alpha$data$sites <- factor(alpha$data$sites, levels = siteorder$`sample_metadata$sites`)
  pdf(file='alphaDiversity_Taxa_normalized_siteClusteredSamples_alltaxa_Chao1Shannon.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
  
  total = median(sample_sums(phylo15b))
  alpha <- plot_richness(transform_sample_counts(phylo15b, standf), nrow = NULL) +
    xlab("Sites")
  alpha$data$samples <- factor(alpha$data$samples, levels = siteorder$`sample_metadata$sites`)
  pdf(file='alphaDiversity_Taxa_normalized_siteGroupedSamples_alltaxa_allMeasures.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
  
  alpha <- plot_richness(transform_sample_counts(phylo15b, standf), measures = c("Chao1","Shannon")) +
    xlab("Sites")
  alpha$data$samples <- factor(alpha$data$samples, levels = siteorder$`sample_metadata$sites`)
  pdf(file='alphaDiversity_Taxa_normalized_siteGroupedSamples_alltaxa_Chao1Shannon.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
}
}, error=function(e){})

tryCatch({
total = median(sample_sums(phylo10))
alpha <- plot_richness(transform_sample_counts(phylo10, standf), nrow = NULL) +
  xlab("Samples")
alpha$data$samples <- factor(alpha$data$samples, levels = sample_order$V1)
pdf(file='alphaDiversity_Taxa_normalized_filtsamples_alltaxa_allMeasures.pdf', width = 15, height = 8.5)
print(alpha)
dev.off()

alpha <- plot_richness(transform_sample_counts(phylo10, standf), measures = c("Chao1","Shannon")) +
  xlab("Samples")
alpha$data$samples <- factor(alpha$data$samples, levels = sample_order$V1)
pdf(file='alphaDiversity_Taxa_normalized_filtsamples_alltaxa_Chao1Shannon.pdf', width = 15, height = 8.5)
print(alpha)
dev.off()
}, error=function(e){})


if (groupingFlag == TRUE) {
  for (i in 1:numberofGroups) {
    tryCatch({
    if (replicateFlag == TRUE) {
      total = median(sample_sums(phylo10))
      alpha <- plot_richness(transform_sample_counts(phylo10, standf), nrow = NULL, x = "replicates", color = paste0("group",i)) +
        xlab("Replicates")
      alpha$data$replicates <- factor(alpha$data$replicates, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_replicateClusteredSamples_alltaxa_allMeasures_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(phylo10, standf), measures = c("Chao1","Shannon"), x = "replicates", color = paste0("group",i)) +
        xlab("Replicates")
      alpha$data$replicates <- factor(alpha$data$replicates, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_replicateClusteredSamples_alltaxa_Chao1Shannon_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      total = median(sample_sums(phylo15a))
      alpha <- plot_richness(transform_sample_counts(phylo15a, standf), nrow = NULL, color = paste0("group",i)) +
        xlab("Replicates")
      alpha$data$samples <- factor(alpha$data$samples, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_replicateGroupedSamples_alltaxa_allMeasures_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(phylo15a, standf), measures = c("Chao1","Shannon"), color = paste0("group",i)) +
        xlab("Replicates")
      alpha$data$samples <- factor(alpha$data$samples, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_replicateGroupedSamples_alltaxa_Chao1Shannon_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    if (sitelabelFlag == TRUE) {
      total = median(sample_sums(phylo10))
      alpha <- plot_richness(transform_sample_counts(phylo10, standf), nrow = NULL, x = "sites", color = paste0("group",i)) +
        xlab("Sites")
      alpha$data$sites <- factor(alpha$data$sites, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_siteClusteredSamples_alltaxa_allMeasures_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(phylo10, standf), measures = c("Chao1","Shannon"), x = "sites", color = paste0("group",i)) +
        xlab("Sites")
      alpha$data$sites <- factor(alpha$data$sites, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_siteClusteredSamples_alltaxa_Chao1Shannon_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      total = median(sample_sums(phylo15b))
      alpha <- plot_richness(transform_sample_counts(phylo15b, standf), nrow = NULL, color = paste0("group",i)) +
        xlab("Sites")
      alpha$data$samples <- factor(alpha$data$samples, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_siteGroupedSamples_alltaxa_allMeasures_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(phylo15b, standf), measures = c("Chao1","Shannon"), color = paste0("group",i)) +
        xlab("Sites")
      alpha$data$samples <- factor(alpha$data$samples, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_siteGroupedSamples_alltaxa_Chao1Shannon_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    total = median(sample_sums(phylo10))
    alpha <- plot_richness(transform_sample_counts(phylo10, standf), nrow = NULL, color = paste0("group",i), sortby = "Chao1") +
      xlab("Samples")
    pdf(file=paste0('alphaDiversity_Taxa_normalized_filtsamples_alltaxa_allMeasures_group',i,'Colored_sortonChao1.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha <- plot_richness(transform_sample_counts(phylo10, standf), nrow = NULL, color = paste0("group",i), sortby = "Shannon") +
      xlab("Samples")
    pdf(file=paste0('alphaDiversity_Taxa_normalized_filtsamples_alltaxa_allMeasures_group',i,'Colored_sortonShannon.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha$data$samples <- factor(alpha$data$samples, levels = sample_order$V1)
    pdf(file=paste0('alphaDiversity_Taxa_normalized_filtsamples_alltaxa_allMeasures_group',i,'Colored.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    
    alpha <- plot_richness(transform_sample_counts(phylo10, standf), nrow = NULL, color = paste0("group",i), measures = c("Chao1","Shannon"), sortby = "Chao1") +
      xlab("Samples")
    pdf(file=paste0('alphaDiversity_Taxa_normalized_filtsamples_alltaxa_Chao1Shannon_group',i,'Colored_sortonChao1.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha <- plot_richness(transform_sample_counts(phylo10, standf), nrow = NULL, color = paste0("group",i), measures = c("Chao1","Shannon"), sortby = "Shannon") +
      xlab("Samples")
    pdf(file=paste0('alphaDiversity_Taxa_normalized_filtsamples_alltaxa_Chao1Shannon_group',i,'Colored_sortonShannon.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha$data$samples <- factor(alpha$data$samples, levels = sample_order$V1)
    pdf(file=paste0('alphaDiversity_Taxa_normalized_filtsamples_alltaxa_Chao1Shannon_group',i,'Colored.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    }, error=function(e){})
  }
}

if (chemDataFlag == TRUE) {
  for (i in chem_headers$V1) {
    tryCatch({
    if (replicateFlag == TRUE) {
      classgroup <- class(sample_data(phylo10)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      total = median(sample_sums(phylo10))
      alpha <- plot_richness(transform_sample_counts(phylo10, standf), nrow = NULL, x = "replicates", color = paste0(i)) +
        xlab("Replicates") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$replicates <- factor(alpha$data$replicates, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_replicateClusteredSamples_alltaxa_allMeasures_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(phylo10, standf), measures = c("Chao1","Shannon"), x = "replicates", color = paste0(i)) +
        xlab("Replicates") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$replicates <- factor(alpha$data$replicates, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_replicateClusteredSamples_alltaxa_Chao1Shannon_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      classgroup <- class(sample_data(phylo15a)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      total = median(sample_sums(phylo15a))
      alpha <- plot_richness(transform_sample_counts(phylo15a, standf), nrow = NULL, color = paste0(i)) +
        xlab("Replicates") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$samples <- factor(alpha$data$samples, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_replicateGroupedSamples_alltaxa_allMeasures_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(phylo15a, standf), measures = c("Chao1","Shannon"), color = paste0(i)) +
        xlab("Replicates") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$samples <- factor(alpha$data$samples, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_replicateGroupedSamples_alltaxa_Chao1Shannon_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    if (sitelabelFlag == TRUE) {
      classgroup <- class(sample_data(phylo10)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      total = median(sample_sums(phylo10))
      alpha <- plot_richness(transform_sample_counts(phylo10, standf), nrow = NULL, x = "sites", color = paste0(i)) +
        xlab("Sites") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$sites <- factor(alpha$data$sites, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_siteClusteredSamples_alltaxa_allMeasures_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(phylo10, standf), measures = c("Chao1","Shannon"), x = "sites", color = paste0(i)) +
        xlab("Sites") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$sites <- factor(alpha$data$sites, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_siteClusteredSamples_alltaxa_Chao1Shannon_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      classgroup <- class(sample_data(phylo15b)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      total = median(sample_sums(phylo15b))
      alpha <- plot_richness(transform_sample_counts(phylo15b, standf), nrow = NULL, color = paste0(i)) +
        xlab("Sites") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$samples <- factor(alpha$data$samples, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_siteGroupedSamples_alltaxa_allMeasures_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(phylo15b, standf), measures = c("Chao1","Shannon"), color = paste0(i)) +
        xlab("Sites") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$samples <- factor(alpha$data$samples, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_siteGroupedSamples_alltaxa_Chao1Shannon_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    classgroup <- class(sample_data(phylo10)[[i]])
    if (classgroup != "integer" && classgroup != "numeric") {
      viridis_scale <- 'scale_color_viridis_d(option="C")'
    } else {
      viridis_scale <- 'scale_color_viridis_c(option="C")'
    }
    
    total = median(sample_sums(phylo10))
    alpha <- plot_richness(transform_sample_counts(phylo10, standf), nrow = NULL, color = paste0(i), sortby = "Chao1") +
      xlab("Samples") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
    pdf(file=paste0('alphaDiversity_Taxa_normalized_filtsamples_alltaxa_allMeasures_chem_',i,'_Colored_sortonChao1.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha <- plot_richness(transform_sample_counts(phylo10, standf), nrow = NULL, color = paste0(i), sortby = "Shannon") +
      xlab("Samples") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
    pdf(file=paste0('alphaDiversity_Taxa_normalized_filtsamples_alltaxa_allMeasures_chem_',i,'_Colored_sortonShannon.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha$data$samples <- factor(alpha$data$samples, levels = sample_order$V1)
    pdf(file=paste0('alphaDiversity_Taxa_normalized_filtsamples_alltaxa_allMeasures_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    
    alpha <- plot_richness(transform_sample_counts(phylo10, standf), nrow = NULL, color = paste0(i), measures = c("Chao1","Shannon"), sortby = "Chao1") +
      xlab("Samples") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
    pdf(file=paste0('alphaDiversity_Taxa_normalized_filtsamples_alltaxa_Chao1Shannon_chem_',i,'_Colored_sortonChao1.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha <- plot_richness(transform_sample_counts(phylo10, standf), nrow = NULL, color = paste0(i), measures = c("Chao1","Shannon"), sortby = "Shannon") +
      xlab("Samples") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
    pdf(file=paste0('alphaDiversity_Taxa_normalized_filtsamples_alltaxa_Chao1Shannon_chem_',i,'_Colored_sortonShannon.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha$data$samples <- factor(alpha$data$samples, levels = sample_order$V1)
    pdf(file=paste0('alphaDiversity_Taxa_normalized_filtsamples_alltaxa_Chao1Shannon_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    }, error=function(e){})
  }
}

###############
# Taxa-based, filter to include only taxa with species level identifiers
###############
tryCatch({
if (replicateFlag == TRUE) {
  total = median(sample_sums(phylo10))
  alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf), nrow = NULL, x = "replicates") +
    xlab("Replicates")
  alpha$data$replicates <- factor(alpha$data$replicates, levels = replicateorder$`sample_metadata$replicates`)
  pdf(file='alphaDiversity_Taxa_normalized_TOSPECIES_replicateClusteredSamples_alltaxa_allMeasures.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
  
  alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf), measures = c("Chao1","Shannon"), x = "replicates") +
    xlab("Replicates")
  alpha$data$replicates <- factor(alpha$data$replicates, levels = replicateorder$`sample_metadata$replicates`)
  pdf(file='alphaDiversity_Taxa_normalized_TOSPECIES_replicateClusteredSamples_alltaxa_Chao1Shannon.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
  
  total = median(sample_sums(phylo15a))
  alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo15a, !is.na(Species)), standf), nrow = NULL) +
    xlab("Replicates")
  alpha$data$samples <- factor(alpha$data$samples, levels = replicateorder$`sample_metadata$replicates`)
  pdf(file='alphaDiversity_Taxa_normalized_TOSPECIES_replicateGroupedSamples_alltaxa_allMeasures.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
  
  alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo15a, !is.na(Species)), standf), measures = c("Chao1","Shannon")) +
    xlab("Replicates")
  alpha$data$samples <- factor(alpha$data$samples, levels = replicateorder$`sample_metadata$replicates`)
  pdf(file='alphaDiversity_Taxa_normalized_TOSPECIES_replicateGroupedSamples_alltaxa_Chao1Shannon.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
}
}, error=function(e){})

tryCatch({
if (sitelabelFlag == TRUE) {
  total = median(sample_sums(phylo10))
  alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf), nrow = NULL, x = "sites") +
    xlab("Sites")
  alpha$data$sites <- factor(alpha$data$sites, levels = siteorder$`sample_metadata$sites`)
  pdf(file='alphaDiversity_Taxa_normalized_TOSPECIES_siteClusteredSamples_alltaxa_allMeasures.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
  
  alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf), measures = c("Chao1","Shannon"), x = "sites") +
    xlab("Sites")
  alpha$data$sites <- factor(alpha$data$sites, levels = siteorder$`sample_metadata$sites`)
  pdf(file='alphaDiversity_Taxa_normalized_TOSPECIES_siteClusteredSamples_alltaxa_Chao1Shannon.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
  
  total = median(sample_sums(phylo15b))
  alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo15b, !is.na(Species)), standf), nrow = NULL) +
    xlab("Sites")
  alpha$data$samples <- factor(alpha$data$samples, levels = siteorder$`sample_metadata$sites`)
  pdf(file='alphaDiversity_Taxa_normalized_TOSPECIES_siteGroupedSamples_alltaxa_allMeasures.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
  
  alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo15b, !is.na(Species)), standf), measures = c("Chao1","Shannon")) +
    xlab("Sites")
  alpha$data$samples <- factor(alpha$data$samples, levels = siteorder$`sample_metadata$sites`)
  pdf(file='alphaDiversity_Taxa_normalized_TOSPECIES_siteGroupedSamples_alltaxa_Chao1Shannon.pdf', width = 15, height = 8.5)
  print(alpha)
  dev.off()
}
}, error=function(e){})

tryCatch({
total = median(sample_sums(phylo10))
alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf), nrow = NULL) +
  xlab("Samples")
alpha$data$samples <- factor(alpha$data$samples, levels = sample_order$V1)
pdf(file='alphaDiversity_Taxa_normalized_TOSPECIES_filtsamples_alltaxa_allMeasures.pdf', width = 15, height = 8.5)
print(alpha)
dev.off()

alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf), measures = c("Chao1","Shannon")) +
  xlab("Samples")
alpha$data$samples <- factor(alpha$data$samples, levels = sample_order$V1)
pdf(file='alphaDiversity_Taxa_normalized_TOSPECIES_filtsamples_alltaxa_Chao1Shannon.pdf', width = 15, height = 8.5)
print(alpha)
dev.off()
}, error=function(e){})


if (groupingFlag == TRUE) {
  for (i in 1:numberofGroups) {
    tryCatch({
    if (replicateFlag == TRUE) {
      total = median(sample_sums(phylo10))
      alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf), nrow = NULL, x = "replicates", color = paste0("group",i)) +
        xlab("Replicates")
      alpha$data$replicates <- factor(alpha$data$replicates, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_replicateClusteredSamples_alltaxa_allMeasures_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf), measures = c("Chao1","Shannon"), x = "replicates", color = paste0("group",i)) +
        xlab("Replicates")
      alpha$data$replicates <- factor(alpha$data$replicates, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_replicateClusteredSamples_alltaxa_Chao1Shannon_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      total = median(sample_sums(phylo15a))
      alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo15a, !is.na(Species)), standf), nrow = NULL, color = paste0("group",i)) +
        xlab("Replicates")
      alpha$data$samples <- factor(alpha$data$samples, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_replicateGroupedSamples_alltaxa_allMeasures_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo15a, !is.na(Species)), standf), measures = c("Chao1","Shannon"), color = paste0("group",i)) +
        xlab("Replicates")
      alpha$data$samples <- factor(alpha$data$samples, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_replicateGroupedSamples_alltaxa_Chao1Shannon_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    if (sitelabelFlag == TRUE) {
      total = median(sample_sums(phylo10))
      alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf), nrow = NULL, x = "sites", color = paste0("group",i)) +
        xlab("Sites")
      alpha$data$sites <- factor(alpha$data$sites, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_siteClusteredSamples_alltaxa_allMeasures_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf), measures = c("Chao1","Shannon"), x = "sites", color = paste0("group",i)) +
        xlab("Sites")
      alpha$data$sites <- factor(alpha$data$sites, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_siteClusteredSamples_alltaxa_Chao1Shannon_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      total = median(sample_sums(phylo15b))
      alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo15b, !is.na(Species)), standf), nrow = NULL, color = paste0("group",i)) +
        xlab("Sites")
      alpha$data$samples <- factor(alpha$data$samples, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_siteGroupedSamples_alltaxa_allMeasures_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo15b, !is.na(Species)), standf), measures = c("Chao1","Shannon"), color = paste0("group",i)) +
        xlab("Sites")
      alpha$data$samples <- factor(alpha$data$samples, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_siteGroupedSamples_alltaxa_Chao1Shannon_group',i,'Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    total = median(sample_sums(phylo10))
    alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf), nrow = NULL, color = paste0("group",i), sortby = "Chao1") +
      xlab("Samples")
    pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_filtsamples_alltaxa_allMeasures_group',i,'Colored_sortonChao1.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf), nrow = NULL, color = paste0("group",i), sortby = "Shannon") +
      xlab("Samples")
    pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_filtsamples_alltaxa_allMeasures_group',i,'Colored_sortonShannon.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha$data$samples <- factor(alpha$data$samples, levels = sample_order$V1)
    pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_filtsamples_alltaxa_allMeasures_group',i,'Colored.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    
    alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf), nrow = NULL, color = paste0("group",i), measures = c("Chao1","Shannon"), sortby = "Chao1") +
      xlab("Samples")
    pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_filtsamples_alltaxa_Chao1Shannon_group',i,'Colored_sortonChao1.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf), nrow = NULL, color = paste0("group",i), measures = c("Chao1","Shannon"), sortby = "Shannon") +
      xlab("Samples")
    pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_filtsamples_alltaxa_Chao1Shannon_group',i,'Colored_sortonShannon.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha$data$samples <- factor(alpha$data$samples, levels = sample_order$V1)
    pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_filtsamples_alltaxa_Chao1Shannon_group',i,'Colored.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    }, error=function(e){})
  }
}

if (chemDataFlag == TRUE) {
  for (i in chem_headers$V1) {
    tryCatch({
    if (replicateFlag == TRUE) {
      classgroup <- class(sample_data(phylo10)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      total = median(sample_sums(phylo10))
      alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf), nrow = NULL, x = "replicates", color = paste0(i)) +
        xlab("Replicates") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$replicates <- factor(alpha$data$replicates, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_replicateClusteredSamples_alltaxa_allMeasures_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf), measures = c("Chao1","Shannon"), x = "replicates", color = paste0(i)) +
        xlab("Replicates") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$replicates <- factor(alpha$data$replicates, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_replicateClusteredSamples_alltaxa_Chao1Shannon_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      classgroup <- class(sample_data(phylo15a)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      total = median(sample_sums(phylo15a))
      alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo15a, !is.na(Species)), standf), nrow = NULL, color = paste0(i)) +
        xlab("Replicates") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$samples <- factor(alpha$data$samples, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_replicateGroupedSamples_alltaxa_allMeasures_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo15a, !is.na(Species)), standf), measures = c("Chao1","Shannon"), color = paste0(i)) +
        xlab("Replicates") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$samples <- factor(alpha$data$samples, levels = replicateorder$`sample_metadata$replicates`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_replicateGroupedSamples_alltaxa_Chao1Shannon_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    if (sitelabelFlag == TRUE) {
      classgroup <- class(sample_data(phylo10)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      total = median(sample_sums(phylo10))
      alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf), nrow = NULL, x = "sites", color = paste0(i)) +
        xlab("Sites") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$sites <- factor(alpha$data$sites, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_siteClusteredSamples_alltaxa_allMeasures_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf), measures = c("Chao1","Shannon"), x = "sites", color = paste0(i)) +
        xlab("Sites") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$sites <- factor(alpha$data$sites, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_siteClusteredSamples_alltaxa_Chao1Shannon_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      classgroup <- class(sample_data(phylo15b)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      total = median(sample_sums(phylo15b))
      alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo15b, !is.na(Species)), standf), nrow = NULL, color = paste0(i)) +
        xlab("Sites") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$samples <- factor(alpha$data$samples, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_siteGroupedSamples_alltaxa_allMeasures_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
      
      alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo15b, !is.na(Species)), standf), measures = c("Chao1","Shannon"), color = paste0(i)) +
        xlab("Sites") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
      alpha$data$samples <- factor(alpha$data$samples, levels = siteorder$`sample_metadata$sites`)
      pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_siteGroupedSamples_alltaxa_Chao1Shannon_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
      print(alpha)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    classgroup <- class(sample_data(phylo10)[[i]])
    if (classgroup != "integer" && classgroup != "numeric") {
      viridis_scale <- 'scale_color_viridis_d(option="C")'
    } else {
      viridis_scale <- 'scale_color_viridis_c(option="C")'
    }
    
    total = median(sample_sums(phylo10))
    alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf), nrow = NULL, color = paste0(i), sortby = "Chao1") +
      xlab("Samples") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
    pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_filtsamples_alltaxa_allMeasures_chem_',i,'_Colored_sortonChao1.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf), nrow = NULL, color = paste0(i), sortby = "Shannon") +
      xlab("Samples") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
    pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_filtsamples_alltaxa_allMeasures_chem_',i,'_Colored_sortonShannon.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha$data$samples <- factor(alpha$data$samples, levels = sample_order$V1)
    pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_filtsamples_alltaxa_allMeasures_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    
    alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf), nrow = NULL, color = paste0(i), measures = c("Chao1","Shannon"), sortby = "Chao1") +
      xlab("Samples") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
    pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_filtsamples_alltaxa_Chao1Shannon_chem_',i,'_Colored_sortonChao1.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha <- plot_richness(transform_sample_counts(subset_taxa(phylo10, !is.na(Species)), standf), nrow = NULL, color = paste0(i), measures = c("Chao1","Shannon"), sortby = "Shannon") +
      xlab("Samples") + eval(parse(text = viridis_scale)) + geom_point(shape=1, color = "black")
    pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_filtsamples_alltaxa_Chao1Shannon_chem_',i,'_Colored_sortonShannon.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    alpha$data$samples <- factor(alpha$data$samples, levels = sample_order$V1)
    pdf(file=paste0('alphaDiversity_Taxa_normalized_TOSPECIES_filtsamples_alltaxa_Chao1Shannon_chem_',i,'_Colored.pdf'), width = 15, height = 8.5)
    print(alpha)
    dev.off()
    }, error=function(e){})
  }
}
}, error=function(e){})
########################
#
# ORDINATION  
#
########################
original_replicateFlag <- FALSE
original_sitelabelFlag <- FALSE
if (replicateFlag == TRUE) {
  test_replicate_column <- sample_metadata$replicates
  test_replicate_column <- unique(test_replicate_column)
  if (length(test_replicate_column) == 1) {
    original_replicateFlag <- replicateFlag
    replicateFlag <- FALSE
  }
}
if (sitelabelFlag == TRUE) {
  test_site_column <- sample_metadata$sites
  test_site_column <- unique(test_site_column)
  if (length(test_site_column) == 1) {
    original_sitelabelFlag <- sitelabelFlag
    sitelabelFlag <- FALSE
  }
}

setwd(paste0(as.character(args[1]), "/05_Ordination/ASV_based/read_count", sep = ""))
options(max.print=1000000)

#Create ordinations for almost all phyloseq objects (NMDS & PCoA)
tryCatch({
message("phylo2 NMDS")
phylo2_NMDS.ord <- ordinate(phylo2, "NMDS")
phylo2_NMDS.ord
phylo2_stress <- phylo2_NMDS.ord$stress
ord_plot <- plot_ordination(phylo2, phylo2_NMDS.ord, type="taxa")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
write.table(ord_tab, "NMDS_ASVbased_taxa_readcount_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
ord_plot <- plot_ordination(phylo2, phylo2_NMDS.ord, type="samples")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "Samples")
ord_tab <- ord_tab %>% select(Samples, NMDS1, NMDS2)
write.table(ord_tab, "NMDS_ASVbased_samples_readcount_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
phylo2_PCoA.ord <- ordinate(phylo2, "PCoA")
message("phylo2 PCoA")
phylo2_PCoA.ord
ord_plot <- plot_ordination(phylo2, phylo2_PCoA.ord, type="taxa")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
write.table(ord_tab, "PCoA_ASVbased_taxa_readcount_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
ord_plot <- plot_ordination(phylo2, phylo2_PCoA.ord, type="samples")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "Samples")
ord_tab <- ord_tab %>% select(Samples, Axis.1, Axis.2)
write.table(ord_tab, "PCoA_ASVbased_samples_readcount_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}, error=function(e){})

tryCatch({
message("phylo3 NMDS")
setwd(paste0(as.character(args[1]), "/05_Ordination/ASV_based/relative_abundance", sep = ""))
phylo3_NMDS.ord <- ordinate(phylo3, "NMDS")
phylo3_NMDS.ord
phylo3_stress <- phylo3_NMDS.ord$stress
ord_plot <- plot_ordination(phylo3, phylo3_NMDS.ord, type="taxa")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
write.table(ord_tab, "NMDS_ASVbased_taxa_relabund_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
ord_plot <- plot_ordination(phylo3, phylo3_NMDS.ord, type="samples")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "Samples")
ord_tab <- ord_tab %>% select(Samples, NMDS1, NMDS2)
write.table(ord_tab, "NMDS_ASVbased_samples_relabund_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
phylo3_PCoA.ord <- ordinate(phylo3, "PCoA")
message("phylo3 PCoA")
phylo3_PCoA.ord
ord_plot <- plot_ordination(phylo3, phylo3_PCoA.ord, type="taxa")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
write.table(ord_tab, "PCoA_ASVbased_taxa_relabund_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
ord_plot <- plot_ordination(phylo3, phylo3_PCoA.ord, type="samples")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "Samples")
ord_tab <- ord_tab %>% select(Samples, Axis.1, Axis.2)
write.table(ord_tab, "PCoA_ASVbased_samples_relabund_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}, error=function(e){})

if (replicateFlag == TRUE) {
  setwd(paste0(as.character(args[1]), "/05_Ordination/ASV_based/relative_abundance", sep = ""))
  
  tryCatch({
  message("phylo7a NMDS")
  phylo7a_NMDS.ord <- ordinate(phylo7a, "NMDS")
  phylo7a_NMDS.ord
  phylo7a_stress <- phylo7a_NMDS.ord$stress
  ord_plot <- plot_ordination(phylo7a, phylo7a_NMDS.ord, type="taxa")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
  write.table(ord_tab, "NMDS_ASVbased_taxa_relabund_replicateGroupedSamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  ord_plot <- plot_ordination(phylo7a, phylo7a_NMDS.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "Replicates")
  ord_tab <- ord_tab %>% select(Replicates, NMDS1, NMDS2)
  write.table(ord_tab, "NMDS_ASVbased_samples_relabund_replicateGroupedSamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  phylo7a_PCoA.ord <- ordinate(phylo7a, "PCoA")
  message("phylo7a PCoA")
  phylo7a_PCoA.ord
  ord_plot <- plot_ordination(phylo7a, phylo7a_PCoA.ord, type="taxa")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
  write.table(ord_tab, "PCoA_ASVbased_taxa_relabund_replicateGroupedSamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  ord_plot <- plot_ordination(phylo7a, phylo7a_PCoA.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "Replicates")
  ord_tab <- ord_tab %>% select(Replicates, Axis.1, Axis.2)
  write.table(ord_tab, "PCoA_ASVbased_samples_relabund_replicateGroupedSamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  }, error=function(e){})
  
  tryCatch({
  message("phylo15a NMDS")
  setwd(paste0(as.character(args[1]), "/05_Ordination/Taxonomy_merge_based/relative_abundance", sep = ""))
  phylo15a_NMDS.ord <- ordinate(phylo15a, "NMDS")
  phylo15a_NMDS.ord
  phylo15a_stress <- phylo15a_NMDS.ord$stress
  ord_plot <- plot_ordination(phylo15a, phylo15a_NMDS.ord, type="taxa")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
  write.table(ord_tab, "NMDS_TAXAbased_taxa_relabund_replicateGroupedSamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  ord_plot <- plot_ordination(phylo15a, phylo15a_NMDS.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "Replicates")
  ord_tab <- ord_tab %>% select(Replicates, NMDS1, NMDS2)
  write.table(ord_tab, "NMDS_TAXAbased_samples_relabund_replicateGroupedSamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  phylo15a_PCoA.ord <- ordinate(phylo15a, "PCoA")
  message("phylo15a PCoA")
  phylo15a_PCoA.ord
  ord_plot <- plot_ordination(phylo15a, phylo15a_PCoA.ord, type="taxa")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
  write.table(ord_tab, "PCoA_TAXAbased_taxa_relabund_replicateGroupedSamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  ord_plot <- plot_ordination(phylo15a, phylo15a_PCoA.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "Replicates")
  ord_tab <- ord_tab %>% select(Replicates, Axis.1, Axis.2)
  write.table(ord_tab, "PCoA_TAXAbased_samples_relabund_replicateGroupedSamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  }, error=function(e){})
}

if (sitelabelFlag == TRUE) {
  setwd(paste0(as.character(args[1]), "/05_Ordination/ASV_based/relative_abundance", sep = ""))
  
  tryCatch({
  message("phylo7b NMDS")
  phylo7b_NMDS.ord <- ordinate(phylo7b, "NMDS")
  phylo7b_NMDS.ord
  phylo7b_stress <- phylo7b_NMDS.ord$stress
  ord_plot <- plot_ordination(phylo7b, phylo7b_NMDS.ord, type="taxa")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
  write.table(ord_tab, "NMDS_ASVbased_taxa_relabund_siteGroupedSamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  ord_plot <- plot_ordination(phylo7b, phylo7b_NMDS.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "Sites")
  ord_tab <- ord_tab %>% select(Sites, NMDS1, NMDS2)
  write.table(ord_tab, "NMDS_ASVbased_samples_relabund_siteGroupedSamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  phylo7b_PCoA.ord <- ordinate(phylo7b, "PCoA")
  message("phylo7b PCoA")
  phylo7b_PCoA.ord
  ord_plot <- plot_ordination(phylo7b, phylo7b_PCoA.ord, type="taxa")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
  write.table(ord_tab, "PCoA_ASVbased_taxa_relabund_siteGroupedSamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  ord_plot <- plot_ordination(phylo7b, phylo7b_PCoA.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "Sites")
  ord_tab <- ord_tab %>% select(Sites, Axis.1, Axis.2)
  write.table(ord_tab, "PCoA_ASVbased_samples_relabund_siteGroupedSamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  }, error=function(e){})
  
  tryCatch({
  message("phylo15b NMDS")
  setwd(paste0(as.character(args[1]), "/05_Ordination/Taxonomy_merge_based/relative_abundance", sep = ""))
  phylo15b_NMDS.ord <- ordinate(phylo15b, "NMDS")
  phylo15b_NMDS.ord
  phylo15b_stress <- phylo15b_NMDS.ord$stress
  ord_plot <- plot_ordination(phylo15b, phylo15b_NMDS.ord, type="taxa")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
  write.table(ord_tab, "NMDS_TAXAbased_taxa_relabund_siteGroupedSamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  ord_plot <- plot_ordination(phylo15b, phylo15b_NMDS.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "Sites")
  ord_tab <- ord_tab %>% select(Sites, NMDS1, NMDS2)
  write.table(ord_tab, "NMDS_TAXAbased_samples_relabund_siteGroupedSamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  phylo15b_PCoA.ord <- ordinate(phylo15b, "PCoA")
  message("phylo15b PCoA")
  phylo15b_PCoA.ord
  ord_plot <- plot_ordination(phylo15b, phylo15b_PCoA.ord, type="taxa")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
  write.table(ord_tab, "PCoA_TAXAbased_taxa_relabund_siteGroupedSamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  ord_plot <- plot_ordination(phylo15b, phylo15b_PCoA.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "Sites")
  ord_tab <- ord_tab %>% select(Sites, Axis.1, Axis.2)
  write.table(ord_tab, "PCoA_TAXAbased_samples_relabund_siteGroupedSamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  }, error=function(e){})
}

tryCatch({
message("phylo10 NMDS")
setwd(paste0(as.character(args[1]), "/05_Ordination/Taxonomy_merge_based/read_count", sep = ""))
phylo10_NMDS.ord <- ordinate(phylo10, "NMDS")
phylo10_NMDS.ord
phylo10_stress <- phylo10_NMDS.ord$stress
ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, type="taxa")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
write.table(ord_tab, "NMDS_TAXAbased_taxa_readcount_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, type="samples")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "Samples")
ord_tab <- ord_tab %>% select(Samples, NMDS1, NMDS2)
write.table(ord_tab, "NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
phylo10_PCoA.ord <- ordinate(phylo10, "PCoA")
message("phylo10 PCoA")
phylo10_PCoA.ord
ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, type="taxa")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
write.table(ord_tab, "PCoA_TAXAbased_taxa_readcount_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, type="samples")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "Samples")
ord_tab <- ord_tab %>% select(Samples, Axis.1, Axis.2)
write.table(ord_tab, "PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}, error=function(e){})

tryCatch({
message("phylo11 NMDS")
setwd(paste0(as.character(args[1]), "/05_Ordination/Taxonomy_merge_based/relative_abundance", sep = ""))
phylo11_NMDS.ord <- ordinate(phylo11, "NMDS")
phylo11_NMDS.ord
phylo11_stress <- phylo11_NMDS.ord$stress
ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, type="taxa")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
write.table(ord_tab, "NMDS_TAXAbased_taxa_relabund_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, type="samples")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "Samples")
ord_tab <- ord_tab %>% select(Samples, NMDS1, NMDS2)
write.table(ord_tab, "NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
phylo11_PCoA.ord <- ordinate(phylo11, "PCoA")
message("phylo11 PCoA")
phylo11_PCoA.ord
ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, type="taxa")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
write.table(ord_tab, "PCoA_TAXAbased_taxa_relabund_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, type="samples")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "Samples")
ord_tab <- ord_tab %>% select(Samples, Axis.1, Axis.2)
write.table(ord_tab, "PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}, error=function(e){})

########################
#
# Convex-hull analysis
#
########################
setwd(paste0(as.character(args[1]), "/05_Ordination/ASV_based/read_count", sep = ""))
if (replicateFlag == TRUE) {
  tryCatch({
  message("REPLICATE Convex hulls - phylo2 - NMDS")
  ord_plot <- plot_ordination(phylo2, phylo2_NMDS.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(replicates, NMDS1, NMDS2)
  ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
  ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
  groupNames <- as.vector(unique(ord_tab$replicates))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
    convex_hull_pts<-chull(keep) #find row numbers of the hull points
    convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
    assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
    CH<-get(paste0("CHxy_",i))
    assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
    print(paste0(i))
    print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
    
    WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
    
    window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
    
    CH<-get(paste0("CHxy_all"))
    assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
    WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
    window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_replicates_NMDS_ASVbased_readcount_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  
  message("REPLICATE Convex hulls - phylo2 - PCoA")
  ord_plot <- plot_ordination(phylo2, phylo2_PCoA.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(replicates, Axis.1, Axis.2)
  ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
  ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
  groupNames <- as.vector(unique(ord_tab$replicates))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      convex_hull_pts<-chull(keep) #find row numbers of the hull points
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
      assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      CH<-get(paste0("CHxy_",i))
      assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      print(paste0(i))
      print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
      
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      
      window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
      
      CH<-get(paste0("CHxy_all"))
      assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_replicates_PCoA_ASVbased_readcount_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  }, error=function(e){})
}

if (sitelabelFlag == TRUE) {
  tryCatch({
  message("SITE Convex hulls - phylo2 - NMDS")
  ord_plot <- plot_ordination(phylo2, phylo2_NMDS.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(sites, NMDS1, NMDS2)
  ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
  ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
  groupNames <- as.vector(unique(ord_tab$sites))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
      convex_hull_pts<-chull(keep) #find row numbers of the hull points
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
      assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
      CH<-get(paste0("CHxy_",i))
      assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      print(paste0(i))
      print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
      
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      
      window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
      
      CH<-get(paste0("CHxy_all"))
      assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_sites_NMDS_ASVbased_readcount_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  
  message("SITE Convex hulls - phylo2 - PCoA")
  ord_plot <- plot_ordination(phylo2, phylo2_PCoA.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(sites, Axis.1, Axis.2)
  ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
  ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
  groupNames <- as.vector(unique(ord_tab$sites))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      convex_hull_pts<-chull(keep) #find row numbers of the hull points
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
      assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      CH<-get(paste0("CHxy_",i))
      assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      print(paste0(i))
      print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
      
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      
      window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
      
      CH<-get(paste0("CHxy_all"))
      assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_sites_PCoA_ASVbased_readcount_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  }, error=function(e){})
}

if (groupingFlag == TRUE) {
  for (i in 1:numberofGroups) {
    tryCatch({
    groupnum <- paste0("group",i)
    classgroup <- class(sample_data(phylo2)[[groupnum]])
    if (classgroup != "integer" && classgroup != "numeric") {
      message(paste0("group",i," Convex hulls - phylo2 - NMDS"))
      ord_plot <- plot_ordination(phylo2, phylo2_NMDS.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0("group",i), NMDS1, NMDS2)
      ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
      ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
      groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        } 
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_NMDS_ASVbased_readcount_filtsamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      
      message(paste0("group",i," Convex hulls - phylo2 - PCoA"))
      ord_plot <- plot_ordination(phylo2, phylo2_PCoA.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0("group",i), Axis.1, Axis.2)
      ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
      ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
      groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        }
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_PCoA_ASVbased_readcount_filtsamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
    }
    }, error=function(e){})
  }
}

if (chemDataFlag == TRUE) {
  for (i in chem_headers$V1) {
    tryCatch({
    classgroup <- class(sample_data(phylo2)[[i]])
    if (classgroup != "integer" && classgroup != "numeric") {
      message(paste0("Chem ",i," Convex hulls - phylo2 - NMDS"))
      ord_plot <- plot_ordination(phylo2, phylo2_NMDS.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0(i), NMDS1, NMDS2)
      ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
      ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
      groupNames <- as.vector(unique(ord_tab[[i]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        } 
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_NMDS_ASVbased_readcount_filtsamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      
      message(paste0("Chem ",i," Convex hulls - phylo2 - PCoA"))
      ord_plot <- plot_ordination(phylo2, phylo2_PCoA.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0(i), Axis.1, Axis.2)
      ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
      ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
      groupNames <- as.vector(unique(ord_tab[[i]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        }
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_PCoA_ASVbased_readcount_filtsamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
    }
    }, error=function(e){})
  }
}

setwd(paste0(as.character(args[1]), "/05_Ordination/ASV_based/relative_abundance", sep = ""))
if (replicateFlag == TRUE) {
  tryCatch({
  message("REPLICATE Convex hulls - phylo3 - NMDS")
  ord_plot <- plot_ordination(phylo3, phylo3_NMDS.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(replicates, NMDS1, NMDS2)
  ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
  ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
  groupNames <- as.vector(unique(ord_tab$replicates))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
      convex_hull_pts<-chull(keep) #find row numbers of the hull points
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
      assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
      CH<-get(paste0("CHxy_",i))
      assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      print(paste0(i))
      print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
      
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      
      window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
      
      CH<-get(paste0("CHxy_all"))
      assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_replicates_NMDS_ASVbased_relabund_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  
  message("REPLICATE Convex hulls - phylo3 - PCoA")
  ord_plot <- plot_ordination(phylo3, phylo3_PCoA.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(replicates, Axis.1, Axis.2)
  ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
  ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
  groupNames <- as.vector(unique(ord_tab$replicates))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      convex_hull_pts<-chull(keep) #find row numbers of the hull points
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
      assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      CH<-get(paste0("CHxy_",i))
      assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      print(paste0(i))
      print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
      
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      
      window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
      
      CH<-get(paste0("CHxy_all"))
      assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_replicates_PCoA_ASVbased_relabund_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  }, error=function(e){})
}

if (sitelabelFlag == TRUE) {
  tryCatch({
  message("SITE Convex hulls - phylo3 - NMDS")
  ord_plot <- plot_ordination(phylo3, phylo3_NMDS.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(sites, NMDS1, NMDS2)
  ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
  ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
  groupNames <- as.vector(unique(ord_tab$sites))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
      convex_hull_pts<-chull(keep) #find row numbers of the hull points
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
      assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
      CH<-get(paste0("CHxy_",i))
      assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      print(paste0(i))
      print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
      
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      
      window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
      
      CH<-get(paste0("CHxy_all"))
      assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_sites_NMDS_ASVbased_relabund_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  
  message("SITE Convex hulls - phylo3 - PCoA")
  ord_plot <- plot_ordination(phylo3, phylo3_PCoA.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(sites, Axis.1, Axis.2)
  ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
  ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
  groupNames <- as.vector(unique(ord_tab$sites))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      convex_hull_pts<-chull(keep) #find row numbers of the hull points
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
      assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      CH<-get(paste0("CHxy_",i))
      assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      print(paste0(i))
      print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
      
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      
      window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
      
      CH<-get(paste0("CHxy_all"))
      assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_sites_PCoA_ASVbased_relabund_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  }, error=function(e){})
}

if (groupingFlag == TRUE) {
  for (i in 1:numberofGroups) {
    tryCatch({
    groupnum <- paste0("group",i)
    classgroup <- class(sample_data(phylo3)[[groupnum]])
    if (classgroup != "integer" && classgroup != "numeric") {
      message(paste0("group",i," Convex hulls - phylo3 - NMDS"))
      ord_plot <- plot_ordination(phylo3, phylo3_NMDS.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0("group",i), NMDS1, NMDS2)
      ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
      ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
      groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        } 
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_NMDS_ASVbased_relabund_filtsamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      
      message(paste0("group",i," Convex hulls - phylo3 - PCoA"))
      ord_plot <- plot_ordination(phylo3, phylo3_PCoA.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0("group",i), Axis.1, Axis.2)
      ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
      ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
      groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        }
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_PCoA_ASVbased_relabund_filtsamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
    }
    }, error=function(e){})
  }
}

if (chemDataFlag == TRUE) {
  for (i in chem_headers$V1) {
    tryCatch({
    classgroup <- class(sample_data(phylo3)[[i]])
    if (classgroup != "integer" && classgroup != "numeric") {
      message(paste0("Chem ",i," Convex hulls - phylo3 - NMDS"))
      ord_plot <- plot_ordination(phylo3, phylo3_NMDS.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0(i), NMDS1, NMDS2)
      ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
      ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
      groupNames <- as.vector(unique(ord_tab[[i]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        } 
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_NMDS_ASVbased_relabund_filtsamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      
      message(paste0("Chem ",i," Convex hulls - phylo3 - PCoA"))
      ord_plot <- plot_ordination(phylo3, phylo3_PCoA.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0(i), Axis.1, Axis.2)
      ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
      ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
      groupNames <- as.vector(unique(ord_tab[[i]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        }
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_PCoA_ASVbased_relabund_filtsamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
    }
    }, error=function(e){})
  }
}

setwd(paste0(as.character(args[1]), "/05_Ordination/Taxonomy_merge_based/read_count", sep = ""))
if (replicateFlag == TRUE) {
  tryCatch({
  message("REPLICATE Convex hulls - phylo10 - NMDS")
  ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(replicates, NMDS1, NMDS2)
  ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
  ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
  groupNames <- as.vector(unique(ord_tab$replicates))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
      convex_hull_pts<-chull(keep) #find row numbers of the hull points
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
      assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
      CH<-get(paste0("CHxy_",i))
      assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      print(paste0(i))
      print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
      
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      
      window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
      
      CH<-get(paste0("CHxy_all"))
      assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_replicates_NMDS_TAXAbased_readcount_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  
  message("REPLICATE Convex hulls - phylo10 - PCoA")
  ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(replicates, Axis.1, Axis.2)
  ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
  ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
  groupNames <- as.vector(unique(ord_tab$replicates))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      convex_hull_pts<-chull(keep) #find row numbers of the hull points
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
      assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      CH<-get(paste0("CHxy_",i))
      assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      print(paste0(i))
      print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
      
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      
      window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
      
      CH<-get(paste0("CHxy_all"))
      assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_replicates_PCoA_TAXAbased_readcount_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  }, error=function(e){})
}

if (sitelabelFlag == TRUE) {
  tryCatch({
  message("SITE Convex hulls - phylo10 - NMDS")
  ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(sites, NMDS1, NMDS2)
  ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
  ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
  groupNames <- as.vector(unique(ord_tab$sites))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
      convex_hull_pts<-chull(keep) #find row numbers of the hull points
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
      assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
      CH<-get(paste0("CHxy_",i))
      assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      print(paste0(i))
      print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
      
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      
      window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
      
      CH<-get(paste0("CHxy_all"))
      assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_sites_NMDS_TAXAbased_readcount_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  
  message("SITE Convex hulls - phylo10 - PCoA")
  ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(sites, Axis.1, Axis.2)
  ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
  ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
  groupNames <- as.vector(unique(ord_tab$sites))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      convex_hull_pts<-chull(keep) #find row numbers of the hull points
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
      assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      CH<-get(paste0("CHxy_",i))
      assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      print(paste0(i))
      print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
      
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      
      window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
      
      CH<-get(paste0("CHxy_all"))
      assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_sites_PCoA_TAXAbased_readcount_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  }, error=function(e){})
}

if (groupingFlag == TRUE) {
  for (i in 1:numberofGroups) {
    tryCatch({
    groupnum <- paste0("group",i)
    classgroup <- class(sample_data(phylo10)[[groupnum]])
    if (classgroup != "integer" && classgroup != "numeric") {
      message(paste0("group",i," Convex hulls - phylo10 - NMDS"))
      ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0("group",i), NMDS1, NMDS2)
      ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
      ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
      groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        } 
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_NMDS_TAXAbased_readcount_filtsamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      
      message(paste0("group",i," Convex hulls - phylo10 - PCoA"))
      ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0("group",i), Axis.1, Axis.2)
      ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
      ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
      groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        }
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_PCoA_TAXAbased_readcount_filtsamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
    }
    }, error=function(e){})
  }
}

if (chemDataFlag == TRUE) {
  for (i in chem_headers$V1) {
    tryCatch({
    classgroup <- class(sample_data(phylo10)[[i]])
    if (classgroup != "integer" && classgroup != "numeric") {
      message(paste0("Chem ",i," Convex hulls - phylo10 - NMDS"))
      ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0(i), NMDS1, NMDS2)
      ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
      ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
      groupNames <- as.vector(unique(ord_tab[[i]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        } 
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_NMDS_TAXAbased_readcount_filtsamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      
      message(paste0("Chem ",i," Convex hulls - phylo10 - PCoA"))
      ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0(i), Axis.1, Axis.2)
      ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
      ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
      groupNames <- as.vector(unique(ord_tab[[i]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        }
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_PCoA_TAXAbased_readcount_filtsamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
    }
    }, error=function(e){})
  }
}

setwd(paste0(as.character(args[1]), "/05_Ordination/Taxonomy_merge_based/relative_abundance", sep = ""))
if (replicateFlag == TRUE) {
  tryCatch({
  message("REPLICATE Convex hulls - phylo11 - NMDS")
  ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(replicates, NMDS1, NMDS2)
  ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
  ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
  groupNames <- as.vector(unique(ord_tab$replicates))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
      convex_hull_pts<-chull(keep) #find row numbers of the hull points
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
      assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
      CH<-get(paste0("CHxy_",i))
      assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      print(paste0(i))
      print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
      
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      
      window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
      
      CH<-get(paste0("CHxy_all"))
      assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_replicates_NMDS_TAXAbased_relabund_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  
  message("REPLICATE Convex hulls - phylo11 - PCoA")
  ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(replicates, Axis.1, Axis.2)
  ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
  ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
  groupNames <- as.vector(unique(ord_tab$replicates))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      convex_hull_pts<-chull(keep) #find row numbers of the hull points
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
      assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      CH<-get(paste0("CHxy_",i))
      assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      print(paste0(i))
      print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
      
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      
      window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
      
      CH<-get(paste0("CHxy_all"))
      assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_replicates_PCoA_TAXAbased_relabund_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  }, error=function(e){})
}

if (sitelabelFlag == TRUE) {
  tryCatch({
  message("SITE Convex hulls - phylo11 - NMDS")
  ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(sites, NMDS1, NMDS2)
  ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
  ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
  groupNames <- as.vector(unique(ord_tab$sites))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
      convex_hull_pts<-chull(keep) #find row numbers of the hull points
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
      assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
      CH<-get(paste0("CHxy_",i))
      assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      print(paste0(i))
      print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
      
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      
      window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
      
      CH<-get(paste0("CHxy_all"))
      assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_sites_NMDS_TAXAbased_relabund_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  
  message("SITE Convex hulls - phylo11 - PCoA")
  ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(sites, Axis.1, Axis.2)
  ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
  ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
  groupNames <- as.vector(unique(ord_tab$sites))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      convex_hull_pts<-chull(keep) #find row numbers of the hull points
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
      assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      CH<-get(paste0("CHxy_",i))
      assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      print(paste0(i))
      print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
      
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      
      window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
      
      CH<-get(paste0("CHxy_all"))
      assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_sites_PCoA_TAXAbased_relabund_filtsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  }, error=function(e){})
}

if (groupingFlag == TRUE) {
  for (i in 1:numberofGroups) {
    tryCatch({
    groupnum <- paste0("group",i)
    classgroup <- class(sample_data(phylo11)[[groupnum]])
    if (classgroup != "integer" && classgroup != "numeric") {
      message(paste0("group",i," Convex hulls - phylo11 - NMDS"))
      ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0("group",i), NMDS1, NMDS2)
      ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
      ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
      groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        } 
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_NMDS_TAXAbased_relabund_filtsamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      
      message(paste0("group",i," Convex hulls - phylo11 - PCoA"))
      ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0("group",i), Axis.1, Axis.2)
      ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
      ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
      groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        }
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_PCoA_TAXAbased_relabund_filtsamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
    }
    }, error=function(e){})
  }
}

if (chemDataFlag == TRUE) {
  for (i in chem_headers$V1) {
    tryCatch({
    classgroup <- class(sample_data(phylo11)[[i]])
    if (classgroup != "integer" && classgroup != "numeric") {
      message(paste0("Chem ",i," Convex hulls - phylo11 - NMDS"))
      ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0(i), NMDS1, NMDS2)
      ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
      ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
      groupNames <- as.vector(unique(ord_tab[[i]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        } 
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_NMDS_TAXAbased_relabund_filtsamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      
      message(paste0("Chem ",i," Convex hulls - phylo11 - PCoA"))
      ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0(i), Axis.1, Axis.2)
      ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
      ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
      groupNames <- as.vector(unique(ord_tab[[i]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        }
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_PCoA_TAXAbased_relabund_filtsamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
    }
    }, error=function(e){})
  }
}

setwd(paste0(as.character(args[1]), "/05_Ordination/ASV_based/relative_abundance", sep = ""))
if (replicateFlag == TRUE) {
  if (groupingFlag == TRUE) {
    for (i in 1:numberofGroups) {
      tryCatch({
      groupnum <- paste0("group",i)
      classgroup <- class(sample_data(phylo7a)[[groupnum]])
      if (classgroup != "integer" && classgroup != "numeric") {
        message(paste0("group",i," Convex hulls - phylo7a - NMDS"))
        ord_plot <- plot_ordination(phylo7a, phylo7a_NMDS.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0("group",i), NMDS1, NMDS2)
        ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
        ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
        groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          } 
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_NMDS_ASVbased_relabund_replicateGroupedSamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
        
        message(paste0("group",i," Convex hulls - phylo7a - PCoA"))
        ord_plot <- plot_ordination(phylo7a, phylo7a_PCoA.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0("group",i), Axis.1, Axis.2)
        ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
        ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
        groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          }
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_PCoA_ASVbased_relabund_replicateGroupedSamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      }
      }, error=function(e){})
    }
  }
  
  if (chemDataFlag == TRUE) {
    for (i in chem_headers$V1) {
      tryCatch({
      classgroup <- class(sample_data(phylo7a)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        message(paste0("Chem ",i," Convex hulls - phylo7a - NMDS"))
        ord_plot <- plot_ordination(phylo7a, phylo7a_NMDS.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0(i), NMDS1, NMDS2)
        ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
        ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
        groupNames <- as.vector(unique(ord_tab[[i]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          } 
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_NMDS_ASVbased_relabund_replicateGroupedSamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
        
        message(paste0("Chem ",i," Convex hulls - phylo7a - PCoA"))
        ord_plot <- plot_ordination(phylo7a, phylo7a_PCoA.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0(i), Axis.1, Axis.2)
        ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
        ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
        groupNames <- as.vector(unique(ord_tab[[i]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          }
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_PCoA_ASVbased_relabund_replicateGroupedSamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      }
      }, error=function(e){})
    }
  }
}

if (sitelabelFlag == TRUE) {
  if (groupingFlag == TRUE) {
    for (i in 1:numberofGroups) {
      tryCatch({
      groupnum <- paste0("group",i)
      classgroup <- class(sample_data(phylo7b)[[groupnum]])
      if (classgroup != "integer" && classgroup != "numeric") {
        message(paste0("group",i," Convex hulls - phylo7b - NMDS"))
        ord_plot <- plot_ordination(phylo7b, phylo7b_NMDS.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0("group",i), NMDS1, NMDS2)
        ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
        ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
        groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          } 
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_NMDS_ASVbased_relabund_siteGroupedSamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
        
        message(paste0("group",i," Convex hulls - phylo7b - PCoA"))
        ord_plot <- plot_ordination(phylo7b, phylo7b_PCoA.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0("group",i), Axis.1, Axis.2)
        ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
        ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
        groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          }
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_PCoA_ASVbased_relabund_siteGroupedSamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      }
      }, error=function(e){})
    }
  }
  
  if (chemDataFlag == TRUE) {
    for (i in chem_headers$V1) {
      tryCatch({
      classgroup <- class(sample_data(phylo7b)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        message(paste0("Chem ",i," Convex hulls - phylo7b - NMDS"))
        ord_plot <- plot_ordination(phylo7b, phylo7b_NMDS.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0(i), NMDS1, NMDS2)
        ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
        ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
        groupNames <- as.vector(unique(ord_tab[[i]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          } 
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_NMDS_ASVbased_relabund_siteGroupedSamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
        
        message(paste0("Chem ",i," Convex hulls - phylo7b - PCoA"))
        ord_plot <- plot_ordination(phylo7b, phylo7b_PCoA.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0(i), Axis.1, Axis.2)
        ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
        ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
        groupNames <- as.vector(unique(ord_tab[[i]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          }
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_PCoA_ASVbased_relabund_siteGroupedSamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      }
      }, error=function(e){})
    }
  }
}

setwd(paste0(as.character(args[1]), "/05_Ordination/Taxonomy_merge_based/relative_abundance", sep = ""))
if (replicateFlag == TRUE) {
  if (groupingFlag == TRUE) {
    for (i in 1:numberofGroups) {
      tryCatch({
      groupnum <- paste0("group",i)
      classgroup <- class(sample_data(phylo15a)[[groupnum]])
      if (classgroup != "integer" && classgroup != "numeric") {
        message(paste0("group",i," Convex hulls - phylo15a - NMDS"))
        ord_plot <- plot_ordination(phylo15a, phylo15a_NMDS.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0("group",i), NMDS1, NMDS2)
        ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
        ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
        groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          } 
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_NMDS_TAXAbased_relabund_replicateGroupedSamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
        
        message(paste0("group",i," Convex hulls - phylo15a - PCoA"))
        ord_plot <- plot_ordination(phylo15a, phylo15a_PCoA.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0("group",i), Axis.1, Axis.2)
        ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
        ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
        groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          }
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_PCoA_TAXAbased_relabund_replicateGroupedSamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      }
      }, error=function(e){})
    }
  }
  
  if (chemDataFlag == TRUE) {
    for (i in chem_headers$V1) {
      tryCatch({
      classgroup <- class(sample_data(phylo15a)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        message(paste0("Chem ",i," Convex hulls - phylo15a - NMDS"))
        ord_plot <- plot_ordination(phylo15a, phylo15a_NMDS.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0(i), NMDS1, NMDS2)
        ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
        ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
        groupNames <- as.vector(unique(ord_tab[[i]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          } 
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_NMDS_TAXAbased_relabund_replicateGroupedSamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
        
        message(paste0("Chem ",i," Convex hulls - phylo15a - PCoA"))
        ord_plot <- plot_ordination(phylo15a, phylo15a_PCoA.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0(i), Axis.1, Axis.2)
        ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
        ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
        groupNames <- as.vector(unique(ord_tab[[i]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          }
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_PCoA_TAXAbased_relabund_replicateGroupedSamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      }
      }, error=function(e){})
    }
  }
}

if (sitelabelFlag == TRUE) {
  if (groupingFlag == TRUE) {
    for (i in 1:numberofGroups) {
      tryCatch({
      groupnum <- paste0("group",i)
      classgroup <- class(sample_data(phylo15b)[[groupnum]])
      if (classgroup != "integer" && classgroup != "numeric") {
        message(paste0("group",i," Convex hulls - phylo15b - NMDS"))
        ord_plot <- plot_ordination(phylo15b, phylo15b_NMDS.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0("group",i), NMDS1, NMDS2)
        ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
        ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
        groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          } 
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_NMDS_TAXAbased_relabund_siteGroupedSamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
        
        message(paste0("group",i," Convex hulls - phylo15b - PCoA"))
        ord_plot <- plot_ordination(phylo15b, phylo15b_PCoA.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0("group",i), Axis.1, Axis.2)
        ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
        ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
        groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          }
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_PCoA_TAXAbased_relabund_siteGroupedSamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      }
      }, error=function(e){})
    }
  }
  
  if (chemDataFlag == TRUE) {
    for (i in chem_headers$V1) {
      tryCatch({
      classgroup <- class(sample_data(phylo15b)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        message(paste0("Chem ",i," Convex hulls - phylo15b - NMDS"))
        ord_plot <- plot_ordination(phylo15b, phylo15b_NMDS.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0(i), NMDS1, NMDS2)
        ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
        ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
        groupNames <- as.vector(unique(ord_tab[[i]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          } 
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_NMDS_TAXAbased_relabund_siteGroupedSamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
        
        message(paste0("Chem ",i," Convex hulls - phylo15b - PCoA"))
        ord_plot <- plot_ordination(phylo15b, phylo15b_PCoA.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0(i), Axis.1, Axis.2)
        ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
        ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
        groupNames <- as.vector(unique(ord_tab[[i]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          }
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_PCoA_TAXAbased_relabund_siteGroupedSamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      }
      }, error=function(e){})
    }
  }
}

########################
#
# ASV-based ORDINATION
#
########################

#ASV-based, read counts (phylo2)
setwd(paste0(as.character(args[1]), "/05_Ordination/ASV_based/read_count", sep = ""))

tryCatch({
ord_plot <- plot_ordination(phylo2, phylo2_NMDS.ord, 
                            type="taxa", 
                            color = "Phylum",
                            title=paste0("NMDS - stress ",round(phylo2_stress, 3)))
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='NMDS_ASVbased_taxa_readcount_filtsamples_alltaxa_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='NMDS_ASVbased_taxa_readcount_filtsamples_alltaxa.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo2, phylo2_NMDS.ord, 
                            type="taxa", 
                            color = "Genus",
                            title=paste0("NMDS - stress ",round(phylo2_stress, 3))) +
  facet_wrap(~Phylum)
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='NMDS_ASVbased_taxa_readcount_filtsamples_alltaxa_facetPhylum_colorGenus_legend.pdf', width = 20, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='NMDS_ASVbased_taxa_readcount_filtsamples_alltaxa_facetPhylum_colorGenus.pdf', width = 15, height = 15)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo2, phylo2_PCoA.ord, 
                            type="taxa", 
                            color = "Phylum",
                            title="PCoA")
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='PCoA_ASVbased_taxa_readcount_filtsamples_alltaxa_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='PCoA_ASVbased_taxa_readcount_filtsamples_alltaxa.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo2, phylo2_PCoA.ord, 
                            type="taxa", 
                            color = "Genus",
                            title="PCoA") +
  facet_wrap(~Phylum)
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='PCoA_ASVbased_taxa_readcount_filtsamples_alltaxa_facetPhylum_colorGenus_legend.pdf', width = 20, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='PCoA_ASVbased_taxa_readcount_filtsamples_alltaxa_facetPhylum_colorGenus.pdf', width = 15, height = 15)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo2, phylo2_NMDS.ord, 
                            type="samples", 
                            title=paste0("NMDS - stress ",round(phylo2_stress, 3))) +
  geom_point(size=3) +
  geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
pdf(file='NMDS_ASVbased_samples_readcount_filtsamples_alltaxa.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo2, phylo2_PCoA.ord, 
                            type="samples", 
                            title="PCoA") +
  geom_point(size=3) +
  geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
pdf(file='PCoA_ASVbased_samples_readcount_filtsamples_alltaxa.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo2, phylo2_NMDS.ord, 
                            type="biplot", 
                            color = "Phylum",
                            title=paste0("NMDS - stress ",round(phylo2_stress, 3)))
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='NMDS_ASVbased_taxa_readcount_filtsamples_alltaxa_biplot_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='NMDS_ASVbased_taxa_readcount_filtsamples_alltaxa_biplot.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo2, phylo2_PCoA.ord, 
                            type="biplot", 
                            color = "Phylum",
                            title="PCoA")
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='PCoA_ASVbased_taxa_readcount_filtsamples_alltaxa_biplot_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='PCoA_ASVbased_taxa_readcount_filtsamples_alltaxa_biplot.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

if (groupingFlag == TRUE) {
  for (i in 1:numberofGroups) {
    tryCatch({
    if (replicateFlag == TRUE) {
      ord_plot <- plot_ordination(phylo2, phylo2_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo2_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('NMDS_ASVbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored_encircleReplicates.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo2, phylo2_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('PCoA_ASVbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored_encircleReplicates.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    if (sitelabelFlag == TRUE) {
      ord_plot <- plot_ordination(phylo2, phylo2_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo2_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group=sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('NMDS_ASVbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored_encircleSites.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo2, phylo2_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group=sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('PCoA_ASVbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored_encircleSites.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    groupnum <- paste0("group",i)
    classgroup <- class(sample_data(phylo2)[[groupnum]])
    if (classgroup == "integer" || classgroup == "numeric") {
      ord_plot <- plot_ordination(phylo2, phylo2_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo2_stress, 3))) +
        geom_point(size=3)
      pdf(file=paste0('NMDS_ASVbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo2, phylo2_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3)
      pdf(file=paste0('PCoA_ASVbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    } else {
      ord_plot <- plot_ordination(phylo2, phylo2_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo2_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group=eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        guides(fill = FALSE)
      pdf(file=paste0('NMDS_ASVbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored_encircleGroup.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo2, phylo2_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group=eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        guides(fill = FALSE)
      pdf(file=paste0('PCoA_ASVbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored_encircleGroup.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
  }
}

if (chemDataFlag == TRUE) {
  for (i in chem_headers$V1) {
    tryCatch({
    if (replicateFlag == TRUE) {
      classgroup <- class(sample_data(phylo2)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      ord_plot <- plot_ordination(phylo2, phylo2_NMDS.ord, 
                                  type="samples",
                                  color=paste0(i),
                                  title=paste0("NMDS - stress ",round(phylo2_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('NMDS_ASVbased_samples_readcount_filtsamples_alltaxa_chem_',i,'_Colored_encircleReplicates.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo2, phylo2_PCoA.ord, 
                                  type="samples", 
                                  color=paste0(i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('PCoA_ASVbased_samples_readcount_filtsamples_alltaxa_chem_',i,'_Colored_encircleReplicates.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    if (sitelabelFlag == TRUE) {
      classgroup <- class(sample_data(phylo2)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      ord_plot <- plot_ordination(phylo2, phylo2_NMDS.ord, 
                                  type="samples",
                                  color=paste0(i),
                                  title=paste0("NMDS - stress ",round(phylo2_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group=sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('NMDS_ASVbased_samples_readcount_filtsamples_alltaxa_chem_',i,'_Colored_encircleSites.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo2, phylo2_PCoA.ord, 
                                  type="samples", 
                                  color=paste0(i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group=sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('PCoA_ASVbased_samples_readcount_filtsamples_alltaxa_chem_',i,'_Colored_encircleSites.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    classgroup <- class(sample_data(phylo2)[[i]])
    if (classgroup != "integer" && classgroup != "numeric") {
      viridis_scale <- 'scale_color_viridis_d(option="C")'
    } else {
      viridis_scale <- 'scale_color_viridis_c(option="C")'
    }
    
    if (classgroup == "integer" || classgroup == "numeric") {
      ord_plot <- plot_ordination(phylo2, phylo2_NMDS.ord, 
                                  type="samples",
                                  color=paste0(i),
                                  title=paste0("NMDS - stress ",round(phylo2_stress, 3))) +
        geom_point(size=3) + eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('NMDS_ASVbased_samples_readcount_filtsamples_alltaxa_chem_',i,'_Colored.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo2, phylo2_PCoA.ord, 
                                  type="samples", 
                                  color=paste0(i),
                                  title="PCoA") +
        geom_point(size=3) + eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('PCoA_ASVbased_samples_readcount_filtsamples_alltaxa_chem_',i,'_Colored.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
  }
}

tryCatch({
if (replicateFlag == TRUE) {
  ord_plot <- plot_ordination(phylo2, phylo2_NMDS.ord, 
                              type="samples",
                              color="replicates",
                              title=paste0("NMDS - stress ",round(phylo2_stress, 3))) +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='NMDS_ASVbased_samples_readcount_filtsamples_alltaxa_encircleReplicates.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  
  ord_plot <- plot_ordination(phylo2, phylo2_PCoA.ord, 
                              type="samples", 
                              color="replicates",
                              title="PCoA") +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='PCoA_ASVbased_samples_readcount_filtsamples_alltaxa_encircleReplicates.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
}
}, error=function(e){})

tryCatch({
if (sitelabelFlag == TRUE) {
  ord_plot <- plot_ordination(phylo2, phylo2_NMDS.ord, 
                              type="samples",
                              color="sites",
                              title=paste0("NMDS - stress ",round(phylo2_stress, 3))) +
    geom_point(size=3) +
    geom_encircle(aes(fill = sites, group=sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='NMDS_ASVbased_samples_readcount_filtsamples_alltaxa_encircleSites.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  
  ord_plot <- plot_ordination(phylo2, phylo2_PCoA.ord, 
                              type="samples", 
                              color="sites",
                              title="PCoA") +
    geom_point(size=3) +
    geom_encircle(aes(fill = sites, group=sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='PCoA_ASVbased_samples_readcount_filtsamples_alltaxa_encircleSites.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
}
}, error=function(e){})

tryCatch({
if (replicateFlag == TRUE && sitelabelFlag == TRUE) {
  ord_plot <- plot_ordination(phylo2, phylo2_NMDS.ord, 
                              type="samples",
                              color="sites",
                              shape="replicates",
                              title=paste0("NMDS - stress ",round(phylo2_stress, 3))) +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='NMDS_ASVbased_samples_readcount_filtsamples_alltaxa_encircleshapeReplicates_colorSites.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  
  ord_plot <- plot_ordination(phylo2, phylo2_PCoA.ord, 
                              type="samples", 
                              color="sites",
                              shape="replicates",
                              title="PCoA") +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='PCoA_ASVbased_samples_readcount_filtsamples_alltaxa_encircleshapeReplicates_colorSites.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
}
}, error=function(e){})

tryCatch({
if (replicateFlag == TRUE) {
  setwd(paste0(as.character(args[1]), "/05_Ordination/ASV_based/relative_abundance", sep = ""))
  
  tryCatch({
  ord_plot <- plot_ordination(phylo7a, phylo7a_NMDS.ord, 
                              type="taxa", 
                              color = "Phylum",
                              title=paste0("NMDS - stress ",round(phylo7a_stress, 3)))
  ord_plot_legend <- get_legend(ord_plot)
  ord_plot <- ord_plot + theme(legend.position='none')
  pdf(file='NMDS_ASVbased_taxa_relabund_replicateGroupedSamples_alltaxa_legend.pdf', width = 11, height = 6)
  print(as_ggplot(ord_plot_legend))
  dev.off()
  pdf(file='NMDS_ASVbased_taxa_relabund_replicateGroupedSamples_alltaxa.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  ord_plot <- plot_ordination(phylo7a, phylo7a_PCoA.ord, 
                              type="taxa", 
                              color = "Phylum",
                              title="PCoA")
  ord_plot_legend <- get_legend(ord_plot)
  ord_plot <- ord_plot + theme(legend.position='none')
  pdf(file='PCoA_ASVbased_taxa_relabund_replicateGroupedSamples_alltaxa_legend.pdf', width = 11, height = 6)
  print(as_ggplot(ord_plot_legend))
  dev.off()
  pdf(file='PCoA_ASVbased_taxa_relabund_replicateGroupedSamples_alltaxa.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  ord_plot <- plot_ordination(phylo7a, phylo7a_NMDS.ord, 
                              type="samples", 
                              title=paste0("NMDS - stress ",round(phylo7a_stress, 3))) +
    geom_point(size=3) +
    geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
  pdf(file='NMDS_ASVbased_samples_relabund_replicateGroupedSamples_alltaxa.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  ord_plot <- plot_ordination(phylo7a, phylo7a_PCoA.ord, 
                              type="samples", 
                              title="PCoA") +
    geom_point(size=3) +
    geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
  pdf(file='PCoA_ASVbased_samples_relabund_replicateGroupedSamples_alltaxa.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  if (groupingFlag == TRUE) {
    for (i in 1:numberofGroups) {
      tryCatch({
      groupnum <- paste0("group",i)
      classgroup <- class(sample_data(phylo7a)[[groupnum]])
      if (classgroup == "integer" || classgroup == "numeric") {
      ord_plot <- plot_ordination(phylo7a, phylo7a_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo7a_stress, 3))) +
        geom_point(size=3)
      pdf(file=paste0('NMDS_ASVbased_samples_relabund_replicateGroupedSamples_alltaxa_group',i,'Colored.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo7a, phylo7a_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3)
      pdf(file=paste0('PCoA_ASVbased_samples_relabund_replicateGroupedSamples_alltaxa_group',i,'Colored.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      } else {
        ord_plot <- plot_ordination(phylo7a, phylo7a_NMDS.ord, 
                                    type="samples",
                                    color=paste0("group",i),
                                    title=paste0("NMDS - stress ",round(phylo7a_stress, 3))) +
          geom_point(size=3) +
          geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group=eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
          guides(fill = FALSE)
        pdf(file=paste0('NMDS_ASVbased_samples_relabund_replicateGroupedSamples_alltaxa_group',i,'Colored_encircleGroup.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
        
        ord_plot <- plot_ordination(phylo7a, phylo7a_PCoA.ord, 
                                    type="samples", 
                                    color=paste0("group",i),
                                    title="PCoA") +
          geom_point(size=3) +
          geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group=eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
          guides(fill = FALSE)
        pdf(file=paste0('PCoA_ASVbased_samples_relabund_replicateGroupedSamples_alltaxa_group',i,'Colored_encircleGroup.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
      }
      }, error=function(e){})
    }
  }
  
  if (chemDataFlag == TRUE) {
    for (i in chem_headers$V1) {
      tryCatch({
      classgroup <- class(sample_data(phylo7a)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      if (classgroup == "integer" || classgroup == "numeric") {
        ord_plot <- plot_ordination(phylo7a, phylo7a_NMDS.ord, 
                                    type="samples",
                                    color=paste0(i),
                                    title=paste0("NMDS - stress ",round(phylo7a_stress, 3))) +
          geom_point(size=3) + eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
        pdf(file=paste0('NMDS_ASVbased_samples_relabund_replicateGroupedSamples_alltaxa_chem_',i,'_Colored.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
        
        ord_plot <- plot_ordination(phylo7a, phylo7a_PCoA.ord, 
                                    type="samples", 
                                    color=paste0(i),
                                    title="PCoA") +
          geom_point(size=3) + eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
        pdf(file=paste0('PCoA_ASVbased_samples_relabund_replicateGroupedSamples_alltaxa_chem_',i,'_Colored.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
      }
      }, error=function(e){})
    }
  }
}
}, error=function(e){})

tryCatch({
if (sitelabelFlag == TRUE) {
  setwd(paste0(as.character(args[1]), "/05_Ordination/ASV_based/relative_abundance", sep = ""))
  
  tryCatch({
  ord_plot <- plot_ordination(phylo7b, phylo7b_NMDS.ord, 
                              type="taxa", 
                              color = "Phylum",
                              title=paste0("NMDS - stress ",round(phylo7b_stress, 3)))
  ord_plot_legend <- get_legend(ord_plot)
  ord_plot <- ord_plot + theme(legend.position='none')
  pdf(file='NMDS_ASVbased_taxa_relabund_siteGroupedSamples_alltaxa_legend.pdf', width = 11, height = 6)
  print(as_ggplot(ord_plot_legend))
  dev.off()
  pdf(file='NMDS_ASVbased_taxa_relabund_siteGroupedSamples_alltaxa.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  ord_plot <- plot_ordination(phylo7b, phylo7b_PCoA.ord, 
                              type="taxa", 
                              color = "Phylum",
                              title="PCoA")
  ord_plot_legend <- get_legend(ord_plot)
  ord_plot <- ord_plot + theme(legend.position='none')
  pdf(file='PCoA_ASVbased_taxa_relabund_siteGroupedSamples_alltaxa_legend.pdf', width = 11, height = 6)
  print(as_ggplot(ord_plot_legend))
  dev.off()
  pdf(file='PCoA_ASVbased_taxa_relabund_siteGroupedSamples_alltaxa.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  ord_plot <- plot_ordination(phylo7b, phylo7b_NMDS.ord, 
                              type="samples", 
                              title=paste0("NMDS - stress ",round(phylo7b_stress, 3))) +
    geom_point(size=3) +
    geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
  pdf(file='NMDS_ASVbased_samples_relabund_siteGroupedSamples_alltaxa.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  ord_plot <- plot_ordination(phylo7b, phylo7b_PCoA.ord, 
                              type="samples", 
                              title="PCoA") +
    geom_point(size=3) +
    geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
  pdf(file='PCoA_ASVbased_samples_relabund_siteGroupedSamples_alltaxa.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  if (groupingFlag == TRUE) {
    for (i in 1:numberofGroups) {
      tryCatch({
      groupnum <- paste0("group",i)
      classgroup <- class(sample_data(phylo7b)[[groupnum]])
      if (classgroup == "integer" || classgroup == "numeric") {
        ord_plot <- plot_ordination(phylo7b, phylo7b_NMDS.ord, 
                                    type="samples",
                                    color=paste0("group",i),
                                    title=paste0("NMDS - stress ",round(phylo7b_stress, 3))) +
          geom_point(size=3)
        pdf(file=paste0('NMDS_ASVbased_samples_relabund_siteGroupedSamples_alltaxa_group',i,'Colored.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
        
        ord_plot <- plot_ordination(phylo7b, phylo7b_PCoA.ord, 
                                    type="samples", 
                                    color=paste0("group",i),
                                    title="PCoA") +
          geom_point(size=3)
        pdf(file=paste0('PCoA_ASVbased_samples_relabund_siteGroupedSamples_alltaxa_group',i,'Colored.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
      } else {
        ord_plot <- plot_ordination(phylo7b, phylo7b_NMDS.ord, 
                                    type="samples",
                                    color=paste0("group",i),
                                    title=paste0("NMDS - stress ",round(phylo7b_stress, 3))) +
          geom_point(size=3) +
          geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group=eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
          guides(fill = FALSE)
        pdf(file=paste0('NMDS_ASVbased_samples_relabund_siteGroupedSamples_alltaxa_group',i,'Colored_encircleGroup.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
        
        ord_plot <- plot_ordination(phylo7b, phylo7b_PCoA.ord, 
                                    type="samples", 
                                    color=paste0("group",i),
                                    title="PCoA") +
          geom_point(size=3) +
          geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group=eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
          guides(fill = FALSE)
        pdf(file=paste0('PCoA_ASVbased_samples_relabund_siteGroupedSamples_alltaxa_group',i,'Colored_encircleGroup.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
      }
      }, error=function(e){})
    }
  }
  
  if (chemDataFlag == TRUE) {
    for (i in chem_headers$V1) {
      tryCatch({
      classgroup <- class(sample_data(phylo7b)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      if (classgroup == "integer" || classgroup == "numeric") {
        ord_plot <- plot_ordination(phylo7b, phylo7b_NMDS.ord, 
                                    type="samples",
                                    color=paste0(i),
                                    title=paste0("NMDS - stress ",round(phylo7b_stress, 3))) +
          geom_point(size=3) + eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
        pdf(file=paste0('NMDS_ASVbased_samples_relabund_siteGroupedSamples_alltaxa_chem_',i,'_Colored.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
        
        ord_plot <- plot_ordination(phylo7b, phylo7b_PCoA.ord, 
                                    type="samples", 
                                    color=paste0(i),
                                    title="PCoA") +
          geom_point(size=3) + eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
        pdf(file=paste0('PCoA_ASVbased_samples_relabund_siteGroupedSamples_alltaxa_chem_',i,'_Colored.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
      }
      }, error=function(e){})
    }
  }
}
}, error=function(e){})

#ASV-based relative abundance (phylo3)
setwd(paste0(as.character(args[1]), "/05_Ordination/ASV_based/relative_abundance", sep = ""))

tryCatch({
ord_plot <- plot_ordination(phylo3, phylo3_NMDS.ord, 
                            type="taxa", 
                            color = "Phylum",
                            title=paste0("NMDS - stress ",round(phylo3_stress, 3)))
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='NMDS_ASVbased_taxa_relabund_filtsamples_alltaxa_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='NMDS_ASVbased_taxa_relabund_filtsamples_alltaxa.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo3, phylo3_NMDS.ord, 
                            type="taxa", 
                            color = "Genus",
                            title=paste0("NMDS - stress ",round(phylo3_stress, 3))) +
  facet_wrap(~Phylum)
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='NMDS_ASVbased_taxa_relabund_filtsamples_alltaxa_facetPhylum_colorGenus_legend.pdf', width = 20, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='NMDS_ASVbased_taxa_relabund_filtsamples_alltaxa_facetPhylum_colorGenus.pdf', width = 15, height = 15)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo3, phylo3_PCoA.ord, 
                            type="taxa", 
                            color = "Phylum",
                            title="PCoA")
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='PCoA_ASVbased_taxa_relabund_filtsamples_alltaxa_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='PCoA_ASVbased_taxa_relabund_filtsamples_alltaxa.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo3, phylo3_PCoA.ord, 
                            type="taxa", 
                            color = "Genus",
                            title="PCoA") +
  facet_wrap(~Phylum)
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='PCoA_ASVbased_taxa_relabund_filtsamples_alltaxa_facetPhylum_colorGenus_legend.pdf', width = 20, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='PCoA_ASVbased_taxa_relabund_filtsamples_alltaxa_facetPhylum_colorGenus.pdf', width = 15, height = 15)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo3, phylo3_NMDS.ord, 
                            type="samples", 
                            title=paste0("NMDS - stress ",round(phylo3_stress, 3))) +
  geom_point(size=3) +
  geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
pdf(file='NMDS_ASVbased_samples_relabund_filtsamples_alltaxa.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo3, phylo3_PCoA.ord, 
                            type="samples", 
                            title="PCoA") +
  geom_point(size=3) +
  geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
pdf(file='PCoA_ASVbased_samples_relabund_filtsamples_alltaxa.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo3, phylo3_NMDS.ord, 
                            type="biplot", 
                            color = "Phylum",
                            title=paste0("NMDS - stress ",round(phylo3_stress, 3)))
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='NMDS_ASVbased_taxa_relabund_filtsamples_alltaxa_biplot_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='NMDS_ASVbased_taxa_relabund_filtsamples_alltaxa_biplot.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo3, phylo3_PCoA.ord, 
                            type="biplot", 
                            color = "Phylum",
                            title="PCoA")
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='PCoA_ASVbased_taxa_relabund_filtsamples_alltaxa_biplot_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='PCoA_ASVbased_taxa_relabund_filtsamples_alltaxa_biplot.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

if (groupingFlag == TRUE) {
  for (i in 1:numberofGroups) {
    tryCatch({
    if (replicateFlag == TRUE) {
      ord_plot <- plot_ordination(phylo3, phylo3_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo3_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('NMDS_ASVbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored_encircleReplicates.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo3, phylo3_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('PCoA_ASVbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored_encircleReplicates.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    if (sitelabelFlag == TRUE) {
      ord_plot <- plot_ordination(phylo3, phylo3_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo3_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group=sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('NMDS_ASVbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored_encircleSites.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo3, phylo3_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group=sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('PCoA_ASVbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored_encircleSites.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    groupnum <- paste0("group",i)
    classgroup <- class(sample_data(phylo3)[[groupnum]])
    if (classgroup == "integer" || classgroup == "numeric") {
      ord_plot <- plot_ordination(phylo3, phylo3_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo3_stress, 3))) +
        geom_point(size=3)
      pdf(file=paste0('NMDS_ASVbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo3, phylo3_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3)
      pdf(file=paste0('PCoA_ASVbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    } else {
      ord_plot <- plot_ordination(phylo3, phylo3_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo3_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group=eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        guides(fill = FALSE)
      pdf(file=paste0('NMDS_ASVbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored_encircleGroup.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo3, phylo3_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group=eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        guides(fill = FALSE)
      pdf(file=paste0('PCoA_ASVbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored_encircleGroup.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
  }
}

if (chemDataFlag == TRUE) {
  for (i in chem_headers$V1) {
    tryCatch({
    if (replicateFlag == TRUE) {
      classgroup <- class(sample_data(phylo3)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      ord_plot <- plot_ordination(phylo3, phylo3_NMDS.ord, 
                                  type="samples",
                                  color=paste0(i),
                                  title=paste0("NMDS - stress ",round(phylo3_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('NMDS_ASVbased_samples_relabund_filtsamples_alltaxa_chem_',i,'_Colored_encircleReplicates.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo3, phylo3_PCoA.ord, 
                                  type="samples", 
                                  color=paste0(i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('PCoA_ASVbased_samples_relabund_filtsamples_alltaxa_chem_',i,'_Colored_encircleReplicates.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    if (sitelabelFlag == TRUE) {
      classgroup <- class(sample_data(phylo3)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      ord_plot <- plot_ordination(phylo3, phylo3_NMDS.ord, 
                                  type="samples",
                                  color=paste0(i),
                                  title=paste0("NMDS - stress ",round(phylo3_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group=sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('NMDS_ASVbased_samples_relabund_filtsamples_alltaxa_chem_',i,'_Colored_encircleSites.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo3, phylo3_PCoA.ord, 
                                  type="samples", 
                                  color=paste0(i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group=sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('PCoA_ASVbased_samples_relabund_filtsamples_alltaxa_chem_',i,'_Colored_encircleSites.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    classgroup <- class(sample_data(phylo3)[[i]])
    if (classgroup != "integer" && classgroup != "numeric") {
      viridis_scale <- 'scale_color_viridis_d(option="C")'
    } else {
      viridis_scale <- 'scale_color_viridis_c(option="C")'
    }
    if (classgroup == "integer" || classgroup == "numeric") {
      ord_plot <- plot_ordination(phylo3, phylo3_NMDS.ord, 
                                  type="samples",
                                  color=paste0(i),
                                  title=paste0("NMDS - stress ",round(phylo3_stress, 3))) +
        geom_point(size=3) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('NMDS_ASVbased_samples_relabund_filtsamples_alltaxa_chem_',i,'_Colored.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo3, phylo3_PCoA.ord, 
                                  type="samples", 
                                  color=paste0(i),
                                  title="PCoA") +
        geom_point(size=3) + eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('PCoA_ASVbased_samples_relabund_filtsamples_alltaxa_chem_',i,'_Colored.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
  }
}

tryCatch({
if (replicateFlag == TRUE) {
  ord_plot <- plot_ordination(phylo3, phylo3_NMDS.ord, 
                              type="samples",
                              color="replicates",
                              title=paste0("NMDS - stress ",round(phylo3_stress, 3))) +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='NMDS_ASVbased_samples_relabund_filtsamples_alltaxa_encircleReplicates.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  
  ord_plot <- plot_ordination(phylo3, phylo3_PCoA.ord, 
                              type="samples", 
                              color="replicates",
                              title="PCoA") +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='PCoA_ASVbased_samples_relabund_filtsamples_alltaxa_encircleReplicates.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
}
}, error=function(e){})

tryCatch({
if (sitelabelFlag == TRUE) {
  ord_plot <- plot_ordination(phylo3, phylo3_NMDS.ord, 
                              type="samples",
                              color="sites",
                              title=paste0("NMDS - stress ",round(phylo3_stress, 3))) +
    geom_point(size=3) +
    geom_encircle(aes(fill = sites, group=sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='NMDS_ASVbased_samples_relabund_filtsamples_alltaxa_encircleSites.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  
  ord_plot <- plot_ordination(phylo3, phylo3_PCoA.ord, 
                              type="samples", 
                              color="sites",
                              title="PCoA") +
    geom_point(size=3) +
    geom_encircle(aes(fill = sites, group=sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='PCoA_ASVbased_samples_relabund_filtsamples_alltaxa_encircleSites.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
}
}, error=function(e){})

tryCatch({
if (replicateFlag == TRUE && sitelabelFlag == TRUE) {
  ord_plot <- plot_ordination(phylo3, phylo3_NMDS.ord, 
                              type="samples",
                              color="sites",
                              shape="replicates",
                              title=paste0("NMDS - stress ",round(phylo3_stress, 3))) +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='NMDS_ASVbased_samples_relabund_filtsamples_alltaxa_encircleshapeReplicates_colorSites.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  
  ord_plot <- plot_ordination(phylo3, phylo3_PCoA.ord, 
                              type="samples", 
                              color="sites",
                              shape="replicates",
                              title="PCoA") +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='PCoA_ASVbased_samples_relabund_filtsamples_alltaxa_encircleshapeReplicates_colorSites.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
}
}, error=function(e){})

########################
#
# Taxa-based ORDINATION
#
########################
#TAXA-based, read counts (phylo10)
setwd(paste0(as.character(args[1]), "/05_Ordination/Taxonomy_merge_based/read_count", sep = ""))

tryCatch({
ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                            type="taxa", 
                            color = "Phylum",
                            title=paste0("NMDS - stress ",round(phylo10_stress, 3)))
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='NMDS_TAXAbased_taxa_readcount_filtsamples_alltaxa_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='NMDS_TAXAbased_taxa_readcount_filtsamples_alltaxa.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                            type="taxa", 
                            color = "Genus",
                            title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
  facet_wrap(~Phylum)
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='NMDS_TAXAbased_taxa_readcount_filtsamples_alltaxa_facetPhylum_colorGenus_legend.pdf', width = 20, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='NMDS_TAXAbased_taxa_readcount_filtsamples_alltaxa_facetPhylum_colorGenus.pdf', width = 15, height = 15)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                            type="taxa", 
                            color = "Phylum",
                            title="PCoA")
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='PCoA_TAXAbased_taxa_readcount_filtsamples_alltaxa_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='PCoA_TAXAbased_taxa_readcount_filtsamples_alltaxa.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                            type="taxa", 
                            color = "Genus",
                            title="PCoA") +
  facet_wrap(~Phylum)
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='PCoA_TAXAbased_taxa_readcount_filtsamples_alltaxa_facetPhylum_colorGenus_legend.pdf', width = 20, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='PCoA_TAXAbased_taxa_readcount_filtsamples_alltaxa_facetPhylum_colorGenus.pdf', width = 15, height = 15)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                            type="samples", 
                            title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
  geom_point(size=3) +
  geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
pdf(file='NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                            type="samples", 
                            title="PCoA") +
  geom_point(size=3) +
  geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
pdf(file='PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                            type="biplot", 
                            color = "Phylum",
                            title=paste0("NMDS - stress ",round(phylo10_stress, 3)))
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='NMDS_TAXAbased_taxa_readcount_filtsamples_alltaxa_biplot_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='NMDS_TAXAbased_taxa_readcount_filtsamples_alltaxa_biplot.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                            type="biplot", 
                            color = "Phylum",
                            title="PCoA")
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='PCoA_TAXAbased_taxa_readcount_filtsamples_alltaxa_biplot_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='PCoA_TAXAbased_taxa_readcount_filtsamples_alltaxa_biplot.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

if (groupingFlag == TRUE) {
  for (i in 1:numberofGroups) {
    tryCatch({
    if (replicateFlag == TRUE) {
      ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored_encircleReplicates.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored_encircleReplicates.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    if (sitelabelFlag == TRUE) {
      ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group=sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored_encircleSites.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group=sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored_encircleSites.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    groupnum <- paste0("group",i)
    classgroup <- class(sample_data(phylo10)[[groupnum]])
    if (classgroup == "integer" || classgroup == "numeric") {
      ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
        geom_point(size=3)
      pdf(file=paste0('NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3)
      pdf(file=paste0('PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    } else {
      ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group = eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        guides(fill = FALSE)
      pdf(file=paste0('NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored_encircleGroup.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group = eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        guides(fill = FALSE)
      pdf(file=paste0('PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored_encircleGroup.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
  }
}

if (chemDataFlag == TRUE) {
  for (i in chem_headers$V1) {
    tryCatch({
    if (replicateFlag == TRUE) {
      classgroup <- class(sample_data(phylo10)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                                  type="samples",
                                  color=paste0(i),
                                  title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa_chem_',i,'_Colored_encircleReplicates.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                                  type="samples", 
                                  color=paste0(i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa_chem_',i,'_Colored_encircleReplicates.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    if (sitelabelFlag == TRUE) {
      classgroup <- class(sample_data(phylo10)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                                  type="samples",
                                  color=paste0(i),
                                  title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group=sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa_chem_',i,'_Colored_encircleSites.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                                  type="samples", 
                                  color=paste0(i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group=sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa_chem_',i,'_Colored_encircleSites.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    classgroup <- class(sample_data(phylo10)[[i]])
    if (classgroup != "integer" && classgroup != "numeric") {
      viridis_scale <- 'scale_color_viridis_d(option="C")'
    } else {
      viridis_scale <- 'scale_color_viridis_c(option="C")'
    }
    
    if (classgroup == "integer" || classgroup == "numeric") {
      ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                                  type="samples",
                                  color=paste0(i),
                                  title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
        geom_point(size=3) + eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa_chem_',i,'_Colored.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                                  type="samples", 
                                  color=paste0(i),
                                  title="PCoA") +
        geom_point(size=3) + eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa_chem_',i,'_Colored.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
  }
}

tryCatch({
if (replicateFlag == TRUE) {
  ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                              type="samples",
                              color="replicates",
                              title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa_encircleReplicates.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  
  ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                              type="samples", 
                              color="replicates",
                              title="PCoA") +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa_encircleReplicates.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
}
}, error=function(e){})

tryCatch({
if (sitelabelFlag == TRUE) {
  ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                              type="samples",
                              color="sites",
                              title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
    geom_point(size=3) +
    geom_encircle(aes(fill = sites, group=sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa_encircleSites.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  
  ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                              type="samples", 
                              color="sites",
                              title="PCoA") +
    geom_point(size=3) +
    geom_encircle(aes(fill = sites, group=sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa_encircleSites.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
}
}, error=function(e){})

tryCatch({
if (replicateFlag == TRUE && sitelabelFlag == TRUE) {
  ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                              type="samples",
                              color="sites",
                              shape="replicates",
                              title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa_encircleshapeReplicates_colorSites.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  
  ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                              type="samples", 
                              color="sites",
                              shape="replicates",
                              title="PCoA") +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa_encircleshapeReplicates_colorSites.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
}
}, error=function(e){})

if (replicateFlag == TRUE) {
  setwd(paste0(as.character(args[1]), "/05_Ordination/Taxonomy_merge_based/relative_abundance", sep = ""))
  
  tryCatch({
  ord_plot <- plot_ordination(phylo15a, phylo15a_NMDS.ord, 
                              type="taxa", 
                              color = "Phylum",
                              title=paste0("NMDS - stress ",round(phylo15a_stress, 3)))
  ord_plot_legend <- get_legend(ord_plot)
  ord_plot <- ord_plot + theme(legend.position='none')
  pdf(file='NMDS_TAXAbased_taxa_relabund_replicateGroupedSamples_alltaxa_legend.pdf', width = 11, height = 6)
  print(as_ggplot(ord_plot_legend))
  dev.off()
  pdf(file='NMDS_TAXAbased_taxa_relabund_replicateGroupedSamples_alltaxa.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  ord_plot <- plot_ordination(phylo15a, phylo15a_PCoA.ord, 
                              type="taxa", 
                              color = "Phylum",
                              title="PCoA")
  ord_plot_legend <- get_legend(ord_plot)
  ord_plot <- ord_plot + theme(legend.position='none')
  pdf(file='PCoA_TAXAbased_taxa_relabund_replicateGroupedSamples_alltaxa_legend.pdf', width = 11, height = 6)
  print(as_ggplot(ord_plot_legend))
  dev.off()
  pdf(file='PCoA_TAXAbased_taxa_relabund_replicateGroupedSamples_alltaxa.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  ord_plot <- plot_ordination(phylo15a, phylo15a_NMDS.ord, 
                              type="samples", 
                              title=paste0("NMDS - stress ",round(phylo15a_stress, 3))) +
    geom_point(size=3) +
    geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
  pdf(file='NMDS_TAXAbased_samples_relabund_replicateGroupedSamples_alltaxa.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  ord_plot <- plot_ordination(phylo15a, phylo15a_PCoA.ord, 
                              type="samples", 
                              title="PCoA") +
    geom_point(size=3) +
    geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
  pdf(file='PCoA_TAXAbased_samples_relabund_replicateGroupedSamples_alltaxa.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  if (groupingFlag == TRUE) {
    for (i in 1:numberofGroups) {
      tryCatch({
      groupnum <- paste0("group",i)
      classgroup <- class(sample_data(phylo15a)[[groupnum]])
      if (classgroup == "integer" || classgroup == "numeric") {
        ord_plot <- plot_ordination(phylo15a, phylo15a_NMDS.ord, 
                                    type="samples",
                                    color=paste0("group",i),
                                    title=paste0("NMDS - stress ",round(phylo15a_stress, 3))) +
          geom_point(size=3)
        pdf(file=paste0('NMDS_TAXAbased_samples_relabund_replicateGroupedSamples_alltaxa_group',i,'Colored.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
        
        ord_plot <- plot_ordination(phylo15a, phylo15a_PCoA.ord, 
                                    type="samples", 
                                    color=paste0("group",i),
                                    title="PCoA") +
          geom_point(size=3)
        pdf(file=paste0('PCoA_TAXAbased_samples_relabund_replicateGroupedSamples_alltaxa_group',i,'Colored.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
      } else {
        ord_plot <- plot_ordination(phylo15a, phylo15a_NMDS.ord, 
                                    type="samples",
                                    color=paste0("group",i),
                                    title=paste0("NMDS - stress ",round(phylo15a_stress, 3))) +
          geom_point(size=3) +
          geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group=eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
          guides(fill = FALSE)
        pdf(file=paste0('NMDS_TAXAbased_samples_relabund_replicateGroupedSamples_alltaxa_group',i,'Colored_encircleGroup.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
        
        ord_plot <- plot_ordination(phylo15a, phylo15a_PCoA.ord, 
                                    type="samples", 
                                    color=paste0("group",i),
                                    title="PCoA") +
          geom_point(size=3) +
          geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group=eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
          guides(fill = FALSE)
        pdf(file=paste0('PCoA_TAXAbased_samples_relabund_replicateGroupedSamples_alltaxa_group',i,'Colored_encircleGroup.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
      }
      }, error=function(e){})
    }
  }
  
  if (chemDataFlag == TRUE) {
    for (i in chem_headers$V1) {
      tryCatch({
      classgroup <- class(sample_data(phylo15a)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      if (classgroup == "integer" || classgroup == "numeric") {
        ord_plot <- plot_ordination(phylo15a, phylo15a_NMDS.ord, 
                                    type="samples",
                                    color=paste0(i),
                                    title=paste0("NMDS - stress ",round(phylo15a_stress, 3))) +
          geom_point(size=3) + eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
        pdf(file=paste0('NMDS_TAXAbased_samples_relabund_replicateGroupedSamples_alltaxa_chem_',i,'_Colored.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
        
        ord_plot <- plot_ordination(phylo15a, phylo15a_PCoA.ord, 
                                    type="samples", 
                                    color=paste0(i),
                                    title="PCoA") +
          geom_point(size=3) + eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
        pdf(file=paste0('PCoA_TAXAbased_samples_relabund_replicateGroupedSamples_alltaxa_chem_',i,'_Colored.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
      }
      }, error=function(e){})
    }
  }
}

if (sitelabelFlag == TRUE) {
  setwd(paste0(as.character(args[1]), "/05_Ordination/Taxonomy_merge_based/relative_abundance", sep = ""))
  
  tryCatch({
  ord_plot <- plot_ordination(phylo15b, phylo15b_NMDS.ord, 
                              type="taxa", 
                              color = "Phylum",
                              title=paste0("NMDS - stress ",round(phylo15b_stress, 3)))
  ord_plot_legend <- get_legend(ord_plot)
  ord_plot <- ord_plot + theme(legend.position='none')
  pdf(file='NMDS_TAXAbased_taxa_relabund_siteGroupedSamples_alltaxa_legend.pdf', width = 11, height = 6)
  print(as_ggplot(ord_plot_legend))
  dev.off()
  pdf(file='NMDS_TAXAbased_taxa_relabund_siteGroupedSamples_alltaxa.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  ord_plot <- plot_ordination(phylo15b, phylo15b_PCoA.ord, 
                              type="taxa", 
                              color = "Phylum",
                              title="PCoA")
  ord_plot_legend <- get_legend(ord_plot)
  ord_plot <- ord_plot + theme(legend.position='none')
  pdf(file='PCoA_TAXAbased_taxa_relabund_siteGroupedSamples_alltaxa_legend.pdf', width = 11, height = 6)
  print(as_ggplot(ord_plot_legend))
  dev.off()
  pdf(file='PCoA_TAXAbased_taxa_relabund_siteGroupedSamples_alltaxa.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  ord_plot <- plot_ordination(phylo15b, phylo15b_NMDS.ord, 
                              type="samples", 
                              title=paste0("NMDS - stress ",round(phylo15b_stress, 3))) +
    geom_point(size=3) +
    geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
  pdf(file='NMDS_TAXAbased_samples_relabund_siteGroupedSamples_alltaxa.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  ord_plot <- plot_ordination(phylo15b, phylo15b_PCoA.ord, 
                              type="samples", 
                              title="PCoA") +
    geom_point(size=3) +
    geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
  pdf(file='PCoA_TAXAbased_samples_relabund_siteGroupedSamples_alltaxa.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  if (groupingFlag == TRUE) {
    for (i in 1:numberofGroups) {
      tryCatch({
      groupnum <- paste0("group",i)
      classgroup <- class(sample_data(phylo15b)[[groupnum]])
      if (classgroup == "integer" || classgroup == "numeric") {
        ord_plot <- plot_ordination(phylo15b, phylo15b_NMDS.ord, 
                                    type="samples",
                                    color=paste0("group",i),
                                    title=paste0("NMDS - stress ",round(phylo15b_stress, 3))) +
          geom_point(size=3)
        pdf(file=paste0('NMDS_TAXAbased_samples_relabund_siteGroupedSamples_alltaxa_group',i,'Colored.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
        
        ord_plot <- plot_ordination(phylo15b, phylo15b_PCoA.ord, 
                                    type="samples", 
                                    color=paste0("group",i),
                                    title="PCoA") +
          geom_point(size=3)
        pdf(file=paste0('PCoA_TAXAbased_samples_relabund_siteGroupedSamples_alltaxa_group',i,'Colored.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
      } else {
        ord_plot <- plot_ordination(phylo15b, phylo15b_NMDS.ord, 
                                    type="samples",
                                    color=paste0("group",i),
                                    title=paste0("NMDS - stress ",round(phylo15b_stress, 3))) +
          geom_point(size=3) +
          geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group=eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
          guides(fill = FALSE)
        pdf(file=paste0('NMDS_TAXAbased_samples_relabund_siteGroupedSamples_alltaxa_group',i,'Colored_encircleGroup.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
        
        ord_plot <- plot_ordination(phylo15b, phylo15b_PCoA.ord, 
                                    type="samples", 
                                    color=paste0("group",i),
                                    title="PCoA") +
          geom_point(size=3) +
          geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group=eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
          guides(fill = FALSE)
        pdf(file=paste0('PCoA_TAXAbased_samples_relabund_siteGroupedSamples_alltaxa_group',i,'Colored_encircleGroup.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
      }
      }, error=function(e){})
    }
  }
  
  if (chemDataFlag == TRUE) {
    for (i in chem_headers$V1) {
      tryCatch({
      classgroup <- class(sample_data(phylo15b)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      if (classgroup == "integer" || classgroup == "numeric") {
        ord_plot <- plot_ordination(phylo15b, phylo15b_NMDS.ord, 
                                    type="samples",
                                    color=paste0(i),
                                    title=paste0("NMDS - stress ",round(phylo15b_stress, 3))) +
          geom_point(size=3) + eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
        pdf(file=paste0('NMDS_TAXAbased_samples_relabund_siteGroupedSamples_alltaxa_chem_',i,'_Colored.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
        
        ord_plot <- plot_ordination(phylo15b, phylo15b_PCoA.ord, 
                                    type="samples", 
                                    color=paste0(i),
                                    title="PCoA") +
          geom_point(size=3) + eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
        pdf(file=paste0('PCoA_TAXAbased_samples_relabund_siteGroupedSamples_alltaxa_chem_',i,'_Colored.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
      }
      }, error=function(e){})
    }
  }
}

#Taxa-based relative abundance (phylo11)
setwd(paste0(as.character(args[1]), "/05_Ordination/Taxonomy_merge_based/relative_abundance", sep = ""))

tryCatch({
ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                            type="taxa", 
                            color = "Phylum",
                            title=paste0("NMDS - stress ",round(phylo11_stress, 3)))
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='NMDS_TAXAbased_taxa_relabund_filtsamples_alltaxa_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='NMDS_TAXAbased_taxa_relabund_filtsamples_alltaxa.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                            type="taxa", 
                            color = "Genus",
                            title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
  facet_wrap(~Phylum)
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='NMDS_TAXAbased_taxa_relabund_filtsamples_alltaxa_facetPhylum_colorGenus_legend.pdf', width = 20, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='NMDS_TAXAbased_taxa_relabund_filtsamples_alltaxa_facetPhylum_colorGenus.pdf', width = 15, height = 15)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                            type="taxa", 
                            color = "Phylum",
                            title="PCoA")
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='PCoA_TAXAbased_taxa_relabund_filtsamples_alltaxa_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='PCoA_TAXAbased_taxa_relabund_filtsamples_alltaxa.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                            type="taxa", 
                            color = "Genus",
                            title="PCoA") +
  facet_wrap(~Phylum)
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='PCoA_TAXAbased_taxa_relabund_filtsamples_alltaxa_facetPhylum_colorGenus_legend.pdf', width = 20, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='PCoA_TAXAbased_taxa_relabund_filtsamples_alltaxa_facetPhylum_colorGenus.pdf', width = 15, height = 15)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                            type="samples", 
                            title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
  geom_point(size=3) +
  geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
pdf(file='NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                            type="samples", 
                            title="PCoA") +
  geom_point(size=3) +
  geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
pdf(file='PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                            type="biplot", 
                            color = "Phylum",
                            title=paste0("NMDS - stress ",round(phylo11_stress, 3)))
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='NMDS_TAXAbased_taxa_relabund_filtsamples_alltaxa_biplot_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='NMDS_TAXAbased_taxa_relabund_filtsamples_alltaxa_biplot.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                            type="biplot", 
                            color = "Phylum",
                            title="PCoA")
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='PCoA_TAXAbased_taxa_relabund_filtsamples_alltaxa_biplot_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='PCoA_TAXAbased_taxa_relabund_filtsamples_alltaxa_biplot.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

if (groupingFlag == TRUE) {
  for (i in 1:numberofGroups) {
    tryCatch({
    if (replicateFlag == TRUE) {
      ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored_encircleReplicates.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored_encircleReplicates.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    if (sitelabelFlag == TRUE) {
      ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group=sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored_encircleSites.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group=sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored_encircleSites.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    groupnum <- paste0("group",i)
    classgroup <- class(sample_data(phylo11)[[groupnum]])
    if (classgroup == "integer" || classgroup == "numeric") {
      ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
        geom_point(size=3)
      pdf(file=paste0('NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3)
      pdf(file=paste0('PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    } else {
      ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group=eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        guides(fill = FALSE)
      pdf(file=paste0('NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored_encircleGroup.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group=eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        guides(fill = FALSE)
      pdf(file=paste0('PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored_encircleGroup.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
  }
}

if (chemDataFlag == TRUE) {
  for (i in chem_headers$V1) {
    tryCatch({
    if (replicateFlag == TRUE) {
      classgroup <- class(sample_data(phylo11)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                                  type="samples",
                                  color=paste0(i),
                                  title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa_chem_',i,'_Colored_encircleReplicates.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                                  type="samples", 
                                  color=paste0(i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa_chem_',i,'_Colored_encircleReplicates.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    if (sitelabelFlag == TRUE) {
      classgroup <- class(sample_data(phylo11)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                                  type="samples",
                                  color=paste0(i),
                                  title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa_chem_',i,'_Colored_encircleSites.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                                  type="samples", 
                                  color=paste0(i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group=sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa_chem_',i,'_Colored_encircleSites.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    classgroup <- class(sample_data(phylo11)[[i]])
    if (classgroup != "integer" && classgroup != "numeric") {
      viridis_scale <- 'scale_color_viridis_d(option="C")'
    } else {
      viridis_scale <- 'scale_color_viridis_c(option="C")'
    }
    if (classgroup == "integer" || classgroup == "numeric") {
      ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                                  type="samples",
                                  color=paste0(i),
                                  title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
        geom_point(size=3) + eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa_chem_',i,'_Colored.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                                  type="samples", 
                                  color=paste0(i),
                                  title="PCoA") +
        geom_point(size=3) + eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa_chem_',i,'_Colored.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
  }
}

tryCatch({
if (replicateFlag == TRUE) {
  ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                              type="samples",
                              color="replicates",
                              title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa_encircleReplicates.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  
  ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                              type="samples", 
                              color="replicates",
                              title="PCoA") +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa_encircleReplicates.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
}
}, error=function(e){})

tryCatch({
if (sitelabelFlag == TRUE) {
  ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                              type="samples",
                              color="sites",
                              title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
    geom_point(size=3) +
    geom_encircle(aes(fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa_encircleSites.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  
  ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                              type="samples", 
                              color="sites",
                              title="PCoA") +
    geom_point(size=3) +
    geom_encircle(aes(fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa_encircleSites.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
}
}, error=function(e){})

tryCatch({
if (replicateFlag == TRUE && sitelabelFlag == TRUE) {
  ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                              type="samples",
                              color="sites",
                              shape="replicates",
                              title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa_encircleshapeReplicates_colorSites.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  
  ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                              type="samples", 
                              color="sites",
                              shape="replicates",
                              title="PCoA") +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa_encircleshapeReplicates_colorSites.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
}
}, error=function(e){})


########################
#
# NETWORK ANALYSIS
#
########################
if (original_replicateFlag == TRUE) {
  replicateFlag <- TRUE
}
if (original_sitelabelFlag == TRUE) {
  sitelabelFlag <- TRUE
}


setwd(paste0(as.character(args[1]), "/06_Network/ASV_based/read_count", sep = ""))

for (i in 1:9) {
  tryCatch({
    taxa_net <- make_network(filter_taxa(phylo2, function(x) sum(x >= 1) > 0, TRUE), type="taxa", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
    net_plot <- plot_network(taxa_net, filter_taxa(phylo2, function(x) sum(x >= 1) > 0, TRUE))
    pdf(file=paste0('Network_ASVbased_taxa_readcount_filtsamples_gte1perctaxa_dist0.',i,'_ASVLabeled.pdf'), width = 8, height = 8)
    print(net_plot)
    dev.off()
  }, error=function(e){})
}


# Taxa example with plot_net (does not graph low distance edges, erroneous)
# net_plot <- plot_net(filter_taxa(phylo3, function(x) sum(x > 1) > 0, TRUE), 
#          type="taxa", 
#          color="Phylum", 
#          point_label="Genus",
#          maxdist = 0.7)
# net_plot_legend <- get_legend(net_plot)
# net_plot <- net_plot + theme(legend.position='none')
# pdf(file='Network_ASVbased_taxa_readcount_filtsamples_alltaxa_dist0.7_incorrectNearEdges_GenusLabeled_legend.pdf', width = 11, height = 6)
# print(as_ggplot(net_plot_legend))
# dev.off()
# pdf(file='Network_ASVbased_taxa_readcount_filtsamples_alltaxa_dist0.7_incorrectNearEdges_GenusLabeled.pdf', width = 8, height = 8)
# print(net_plot)
# dev.off()

for (i in 1:9) {
  tryCatch({
    sample_net <- make_network(phylo2, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
    net_plot <- plot_network(sample_net, phylo2)
    pdf(file=paste0('Network_ASVbased_samples_readcount_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled.pdf'), width = 8, height = 8)
    print(net_plot)
    dev.off()
  }, error=function(e){})
}

if (replicateFlag == TRUE) {
  for (i in 1:9) {
    tryCatch({
      sample_net <- make_network(phylo2, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
      net_plot <- plot_network(sample_net, phylo2, color="replicates")
      pdf(file=paste0('Network_ASVbased_samples_readcount_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled_colorReplicates.pdf'), width = 8, height = 8)
      print(net_plot)
      dev.off()
    }, error=function(e){})
  }
}

if (sitelabelFlag == TRUE) {
  for (i in 1:9) {
    tryCatch({
      sample_net <- make_network(phylo2, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
      net_plot <- plot_network(sample_net, phylo2, color="sites")
      pdf(file=paste0('Network_ASVbased_samples_readcount_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled_colorSites.pdf'), width = 8, height = 8)
      print(net_plot)
      dev.off()
    }, error=function(e){})
  }
}

# Sample example with plot_net (does not graph low distance edges, erroneous)
# plot_net(phylo3,
#          type="samples",
#          color="replicates",
#          maxdist = 0.4)


setwd(paste0(as.character(args[1]), "/06_Network/ASV_based/relative_abundance", sep = ""))

for (i in 1:9) {
  tryCatch({
    taxa_net <- make_network(filter_taxa(phylo3, function(x) sum(x >= 1) > 0, TRUE), type="taxa", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
    net_plot <- plot_network(taxa_net, filter_taxa(phylo3, function(x) sum(x >= 1) > 0, TRUE))
    pdf(file=paste0('Network_ASVbased_taxa_relabund_filtsamples_gte1perctaxa_dist0.',i,'_ASVLabeled.pdf'), width = 8, height = 8)
    print(net_plot)
    dev.off()
  }, error=function(e){})
}

for (i in 1:9) {
  tryCatch({
    sample_net <- make_network(phylo3, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
    net_plot <- plot_network(sample_net, phylo3)
    pdf(file=paste0('Network_ASVbased_samples_relabund_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled.pdf'), width = 8, height = 8)
    print(net_plot)
    dev.off()
  }, error=function(e){})
}

if (replicateFlag == TRUE) {
  for (i in 1:9) {
    tryCatch({
      sample_net <- make_network(phylo3, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
      net_plot <- plot_network(sample_net, phylo3, color="replicates")
      pdf(file=paste0('Network_ASVbased_samples_relabund_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled_colorReplicates.pdf'), width = 8, height = 8)
      print(net_plot)
      dev.off()
    }, error=function(e){})
  }
}

if (sitelabelFlag == TRUE) {
  for (i in 1:9) {
    tryCatch({
      sample_net <- make_network(phylo3, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
      net_plot <- plot_network(sample_net, phylo3, color="sites")
      pdf(file=paste0('Network_ASVbased_samples_relabund_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled_colorSites.pdf'), width = 8, height = 8)
      print(net_plot)
      dev.off()
    }, error=function(e){})
  }
}


setwd(paste0(as.character(args[1]), "/06_Network/Taxonomy_merge_based/read_count", sep = ""))

for (i in 1:9) {
  tryCatch({
    taxa_net <- make_network(filter_taxa(phylo10, function(x) sum(x >= 1) > 0, TRUE), type="taxa", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
    net_plot <- plot_network(taxa_net, filter_taxa(phylo10, function(x) sum(x >= 1) > 0, TRUE))
    pdf(file=paste0('Network_TAXAbased_taxa_readcount_filtsamples_gte1perctaxa_dist0.',i,'_ASVLabeled.pdf'), width = 8, height = 8)
    print(net_plot)
    dev.off()
  }, error=function(e){})
}

for (i in 1:9) {
  tryCatch({
    sample_net <- make_network(phylo10, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
    net_plot <- plot_network(sample_net, phylo10)
    pdf(file=paste0('Network_TAXAbased_samples_readcount_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled.pdf'), width = 8, height = 8)
    print(net_plot)
    dev.off()
  }, error=function(e){})
}

if (replicateFlag == TRUE) {
  for (i in 1:9) {
    tryCatch({
      sample_net <- make_network(phylo10, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
      net_plot <- plot_network(sample_net, phylo10, color="replicates")
      pdf(file=paste0('Network_TAXAbased_samples_readcount_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled_colorReplicates.pdf'), width = 8, height = 8)
      print(net_plot)
      dev.off()
    }, error=function(e){})
  }
}

if (sitelabelFlag == TRUE) {
  for (i in 1:9) {
    tryCatch({
      sample_net <- make_network(phylo10, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
      net_plot <- plot_network(sample_net, phylo10, color="sites")
      pdf(file=paste0('Network_TAXAbased_samples_readcount_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled_colorSites.pdf'), width = 8, height = 8)
      print(net_plot)
      dev.off()
    }, error=function(e){})
  }
}


setwd(paste0(as.character(args[1]), "/06_Network/Taxonomy_merge_based/relative_abundance", sep = ""))

for (i in 1:9) {
  tryCatch({
    taxa_net <- make_network(filter_taxa(phylo11, function(x) sum(x >= 1) > 0, TRUE), type="taxa", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
    net_plot <- plot_network(taxa_net, filter_taxa(phylo11, function(x) sum(x >= 1) > 0, TRUE))
    pdf(file=paste0('Network_TAXAbased_taxa_relabund_filtsamples_gte1perctaxa_dist0.',i,'_ASVLabeled.pdf'), width = 8, height = 8)
    print(net_plot)
    dev.off()
  }, error=function(e){})
}

for (i in 1:9) {
  tryCatch({
    sample_net <- make_network(phylo11, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
    net_plot <- plot_network(sample_net, phylo11)
    pdf(file=paste0('Network_TAXAbased_samples_relabund_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled.pdf'), width = 8, height = 8)
    print(net_plot)
    dev.off()
  }, error=function(e){})
}

if (replicateFlag == TRUE) {
  for (i in 1:9) {
    tryCatch({
      sample_net <- make_network(phylo11, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
      net_plot <- plot_network(sample_net, phylo11, color="replicates")
      pdf(file=paste0('Network_TAXAbased_samples_relabund_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled_colorReplicates.pdf'), width = 8, height = 8)
      print(net_plot)
      dev.off()
    }, error=function(e){})
  }
}

if (sitelabelFlag == TRUE) {
  for (i in 1:9) {
    tryCatch({
      sample_net <- make_network(phylo11, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
      net_plot <- plot_network(sample_net, phylo11, color="sites")
      pdf(file=paste0('Network_TAXAbased_samples_relabund_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled_colorSites.pdf'), width = 8, height = 8)
      print(net_plot)
      dev.off()
    }, error=function(e){})
  }
}


########################
#
# Rarefaction curves
#
########################
setwd(paste0(as.character(args[1]), "/07_Rarefaction_Curves", sep = ""))

tryCatch({
pdf(file='RarefactionCurve_ASVbased_readcount_allsamples_alltaxa.pdf', width = 11, height = 8.5)
rarecurve(t(otu_table(phylo1)), step = 50, cex = 0.5)
dev.off()
}, error=function(e){})

tryCatch({
pdf(file='RarefactionCurve_ASVbased_readcount_filtsamples_alltaxa.pdf', width = 11, height = 8.5)
rarecurve(t(otu_table(phylo2)), step = 50, cex = 0.5)
dev.off()
}, error=function(e){})

tryCatch({
pdf(file='RarefactionCurve_TAXAbased_readcount_filtsamples_alltaxa.pdf', width = 11, height = 8.5)
rarecurve(t(otu_table(phylo10)), step = 50, cex = 0.5)
dev.off()
}, error=function(e){})





########################################################################################################################
#
#
# WARNING: EVERYTHING PAST HERE REQUIRES MODIFICATIONS OF ORIGINAL phylo OBJECTS
#
#
########################################################################################################################













if (original_replicateFlag == TRUE) {
  replicateFlag <- FALSE
}
if (original_sitelabelFlag == TRUE) {
  sitelabelFlag <- FALSE
}


########################  
#
# Taxa-based ORDINATION, filtering to include only SPECIES LEVEL IDENTIFICATION
#
########################
#Create ordinations for new filtered phyloseq objects (NMDS & PCoA)

tryCatch({
if (replicateFlag == TRUE) {
  setwd(paste0(as.character(args[1]), "/05_Ordination/Taxonomy_merge_based/filterInclude_TOSPECIES_only/relative_abundance", sep = ""))
  
  message("phylo15a NMDS -- SPECIES filtered")
  phylo15a <- subset_taxa(phylo15a, !is.na(Species))
  phylo15a_NMDS.ord <- ordinate(phylo15a, "NMDS")
  phylo15a_NMDS.ord
  phylo15a_stress <- phylo15a_NMDS.ord$stress
  ord_plot <- plot_ordination(phylo15a, phylo15a_NMDS.ord, type="taxa")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
  write.table(ord_tab, "NMDS_TAXAbased_taxa_relabund_replicateGroupedSamples_alltaxa_filterTOSPECIES.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  ord_plot <- plot_ordination(phylo15a, phylo15a_NMDS.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "Replicates")
  ord_tab <- ord_tab %>% select(Replicates, NMDS1, NMDS2)
  write.table(ord_tab, "NMDS_TAXAbased_samples_relabund_replicateGroupedSamples_alltaxa_filterTOSPECIES.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  phylo15a_PCoA.ord <- ordinate(phylo15a, "PCoA")
  message("phylo15a PCoA -- SPECIES filtered")
  phylo15a_PCoA.ord
  ord_plot <- plot_ordination(phylo15a, phylo15a_PCoA.ord, type="taxa")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
  write.table(ord_tab, "PCoA_TAXAbased_taxa_relabund_replicateGroupedSamples_alltaxa_filterTOSPECIES.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  ord_plot <- plot_ordination(phylo15a, phylo15a_PCoA.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "Replicates")
  ord_tab <- ord_tab %>% select(Replicates, Axis.1, Axis.2)
  write.table(ord_tab, "PCoA_TAXAbased_samples_relabund_replicateGroupedSamples_alltaxa_filterTOSPECIES.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}
}, error=function(e){})

tryCatch({
if (sitelabelFlag == TRUE) {
  setwd(paste0(as.character(args[1]), "/05_Ordination/Taxonomy_merge_based/filterInclude_TOSPECIES_only/relative_abundance", sep = ""))
  
  message("phylo15b NMDS -- SPECIES filtered")
  phylo15b <- subset_taxa(phylo15b, !is.na(Species))
  phylo15b_NMDS.ord <- ordinate(phylo15b, "NMDS")
  phylo15b_NMDS.ord
  phylo15b_stress <- phylo15b_NMDS.ord$stress
  ord_plot <- plot_ordination(phylo15b, phylo15b_NMDS.ord, type="taxa")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
  write.table(ord_tab, "NMDS_TAXAbased_taxa_relabund_siteGroupedSamples_alltaxa_filterTOSPECIES.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  ord_plot <- plot_ordination(phylo15b, phylo15b_NMDS.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "Sites")
  ord_tab <- ord_tab %>% select(Sites, NMDS1, NMDS2)
  write.table(ord_tab, "NMDS_TAXAbased_samples_relabund_siteGroupedSamples_alltaxa_filterTOSPECIES.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  phylo15b_PCoA.ord <- ordinate(phylo15b, "PCoA")
  message("phylo15b PCoA -- SPECIES filtered")
  phylo15b_PCoA.ord
  ord_plot <- plot_ordination(phylo15b, phylo15b_PCoA.ord, type="taxa")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
  write.table(ord_tab, "PCoA_TAXAbased_taxa_relabund_siteGroupedSamples_alltaxa_filterTOSPECIES.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  ord_plot <- plot_ordination(phylo15b, phylo15b_PCoA.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- tibble::rownames_to_column(ord_tab, "Sites")
  ord_tab <- ord_tab %>% select(Sites, Axis.1, Axis.2)
  write.table(ord_tab, "PCoA_TAXAbased_samples_relabund_siteGroupedSamples_alltaxa_filterTOSPECIES.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}
}, error=function(e){})

tryCatch({
message("phylo10 NMDS -- SPECIES filtered")
setwd(paste0(as.character(args[1]), "/05_Ordination/Taxonomy_merge_based/filterInclude_TOSPECIES_only/read_count", sep = ""))
phylo10 <- subset_taxa(phylo10, !is.na(Species))
phylo10_NMDS.ord <- ordinate(phylo10, "NMDS")
phylo10_NMDS.ord
phylo10_stress <- phylo10_NMDS.ord$stress
ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, type="taxa")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
write.table(ord_tab, "NMDS_TAXAbased_taxa_readcount_filtsamples_alltaxa_filterTOSPECIES.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, type="samples")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "Samples")
ord_tab <- ord_tab %>% select(Samples, NMDS1, NMDS2)
write.table(ord_tab, "NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa_filterTOSPECIES.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
phylo10_PCoA.ord <- ordinate(phylo10, "PCoA")
message("phylo10 PCoA -- SPECIES filtered")
phylo10_PCoA.ord
ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, type="taxa")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
write.table(ord_tab, "PCoA_TAXAbased_taxa_readcount_filtsamples_alltaxa_filterTOSPECIES.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, type="samples")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "Samples")
ord_tab <- ord_tab %>% select(Samples, Axis.1, Axis.2)
write.table(ord_tab, "PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa_filterTOSPECIES.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}, error=function(e){})

tryCatch({
message("phylo11 NMDS -- SPECIES filtered")
setwd(paste0(as.character(args[1]), "/05_Ordination/Taxonomy_merge_based/filterInclude_TOSPECIES_only/relative_abundance", sep = ""))
phylo11 <- subset_taxa(phylo11, !is.na(Species))
phylo11_NMDS.ord <- ordinate(phylo11, "NMDS")
phylo11_NMDS.ord
phylo11_stress <- phylo11_NMDS.ord$stress
ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, type="taxa")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
write.table(ord_tab, "NMDS_TAXAbased_taxa_relabund_filtsamples_alltaxa_filterTOSPECIES.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, type="samples")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "Samples")
ord_tab <- ord_tab %>% select(Samples, NMDS1, NMDS2)
write.table(ord_tab, "NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa_filterTOSPECIES.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
phylo11_PCoA.ord <- ordinate(phylo11, "PCoA")
message("phylo11 PCoA -- SPECIES filtered")
phylo11_PCoA.ord
ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, type="taxa")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
write.table(ord_tab, "PCoA_TAXAbased_taxa_relabund_filtsamples_alltaxa_filterTOSPECIES.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, type="samples")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "Samples")
ord_tab <- ord_tab %>% select(Samples, Axis.1, Axis.2)
write.table(ord_tab, "PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa_filterTOSPECIES.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}, error=function(e){})

########################
#
# Convex-hull analysis (filtered phylo objects)
#
########################

setwd(paste0(as.character(args[1]), "/05_Ordination/Taxonomy_merge_based/filterInclude_TOSPECIES_only/read_count", sep = ""))
if (replicateFlag == TRUE) {
  tryCatch({
  message("REPLICATE Convex hulls - phylo10 - NMDS")
  ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(replicates, NMDS1, NMDS2)
  ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
  ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
  groupNames <- as.vector(unique(ord_tab$replicates))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
      convex_hull_pts<-chull(keep) #find row numbers of the hull points
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
      assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
      CH<-get(paste0("CHxy_",i))
      assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      print(paste0(i))
      print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
      
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      
      window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
      
      CH<-get(paste0("CHxy_all"))
      assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_replicates_NMDS_TAXAbased_readcount_filtsamples_filterIncludeToSPECIESONLY.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  
  message("REPLICATE Convex hulls - phylo10 - PCoA")
  ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(replicates, Axis.1, Axis.2)
  ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
  ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
  groupNames <- as.vector(unique(ord_tab$replicates))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      convex_hull_pts<-chull(keep) #find row numbers of the hull points
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
      assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      CH<-get(paste0("CHxy_",i))
      assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      print(paste0(i))
      print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
      
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      
      window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
      
      CH<-get(paste0("CHxy_all"))
      assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_replicates_PCoA_TAXAbased_readcount_filtsamples_filterIncludeToSPECIESONLY.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  }, error=function(e){})
}

if (sitelabelFlag == TRUE) {
  tryCatch({
  message("SITE Convex hulls - phylo10 - NMDS")
  ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(sites, NMDS1, NMDS2)
  ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
  ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
  groupNames <- as.vector(unique(ord_tab$sites))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
      convex_hull_pts<-chull(keep) #find row numbers of the hull points
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
      assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
      CH<-get(paste0("CHxy_",i))
      assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      print(paste0(i))
      print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
      
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      
      window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
      
      CH<-get(paste0("CHxy_all"))
      assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_sites_NMDS_TAXAbased_readcount_filtsamples_filterIncludeToSPECIESONLY.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  
  message("SITE Convex hulls - phylo10 - PCoA")
  ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(sites, Axis.1, Axis.2)
  ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
  ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
  groupNames <- as.vector(unique(ord_tab$sites))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      convex_hull_pts<-chull(keep) #find row numbers of the hull points
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
      assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      CH<-get(paste0("CHxy_",i))
      assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      print(paste0(i))
      print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
      
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      
      window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
      
      CH<-get(paste0("CHxy_all"))
      assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_sites_PCoA_TAXAbased_readcount_filtsamples_filterIncludeToSPECIESONLY.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  }, error=function(e){})
}

if (groupingFlag == TRUE) {
  for (i in 1:numberofGroups) {
    tryCatch({
    groupnum <- paste0("group",i)
    classgroup <- class(sample_data(phylo10)[[groupnum]])
    if (classgroup != "integer" && classgroup != "numeric") {
      message(paste0("group",i," Convex hulls - phylo10 - NMDS"))
      ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0("group",i), NMDS1, NMDS2)
      ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
      ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
      groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        } 
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_NMDS_TAXAbased_readcount_filtsamples_filterIncludeToSPECIESONLY.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      
      message(paste0("group",i," Convex hulls - phylo10 - PCoA"))
      ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0("group",i), Axis.1, Axis.2)
      ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
      ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
      groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        }
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_PCoA_TAXAbased_readcount_filtsamples_filterIncludeToSPECIESONLY.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
    }
    }, error=function(e){})
  }
}

if (chemDataFlag == TRUE) {
  for (i in chem_headers$V1) {
    tryCatch({
    classgroup <- class(sample_data(phylo10)[[i]])
    if (classgroup != "integer" && classgroup != "numeric") {
      message(paste0("Chem ",i," Convex hulls - phylo10 - NMDS"))
      ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0(i), NMDS1, NMDS2)
      ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
      ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
      groupNames <- as.vector(unique(ord_tab[[i]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        } 
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_NMDS_TAXAbased_readcount_filtsamples_filterIncludeToSPECIESONLY.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      
      message(paste0("Chem ",i," Convex hulls - phylo10 - PCoA"))
      ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0(i), Axis.1, Axis.2)
      ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
      ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
      groupNames <- as.vector(unique(ord_tab[[i]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        }
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_PCoA_TAXAbased_readcount_filtsamples_filterIncludeToSPECIESONLY.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
    }
    }, error=function(e){})
  }
}


setwd(paste0(as.character(args[1]), "/05_Ordination/Taxonomy_merge_based/filterInclude_TOSPECIES_only/relative_abundance", sep = ""))
if (replicateFlag == TRUE) {
  tryCatch({
  message("REPLICATE Convex hulls - phylo11 - NMDS")
  ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(replicates, NMDS1, NMDS2)
  ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
  ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
  groupNames <- as.vector(unique(ord_tab$replicates))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
      convex_hull_pts<-chull(keep) #find row numbers of the hull points
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
      assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
      CH<-get(paste0("CHxy_",i))
      assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      print(paste0(i))
      print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
      
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      
      window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
      
      CH<-get(paste0("CHxy_all"))
      assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_replicates_NMDS_TAXAbased_relabund_filtsamples_filterIncludeToSPECIESONLY.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  
  message("REPLICATE Convex hulls - phylo11 - PCoA")
  ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(replicates, Axis.1, Axis.2)
  ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
  ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
  groupNames <- as.vector(unique(ord_tab$replicates))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      convex_hull_pts<-chull(keep) #find row numbers of the hull points
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
      assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      CH<-get(paste0("CHxy_",i))
      assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      print(paste0(i))
      print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
      
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      
      window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
      
      CH<-get(paste0("CHxy_all"))
      assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_replicates_PCoA_TAXAbased_relabund_filtsamples_filterIncludeToSPECIESONLY.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  }, error=function(e){})
}

if (sitelabelFlag == TRUE) {
  tryCatch({
  message("SITE Convex hulls - phylo11 - NMDS")
  ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(sites, NMDS1, NMDS2)
  ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
  ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
  groupNames <- as.vector(unique(ord_tab$sites))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
      convex_hull_pts<-chull(keep) #find row numbers of the hull points
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
      assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$NMDS1) > 2) {
      CH<-get(paste0("CHxy_",i))
      assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      print(paste0(i))
      print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
      
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      
      window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
      
      CH<-get(paste0("CHxy_all"))
      assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_sites_NMDS_TAXAbased_relabund_filtsamples_filterIncludeToSPECIESONLY.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  
  message("SITE Convex hulls - phylo11 - PCoA")
  ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(sites, Axis.1, Axis.2)
  ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
  ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
  groupNames <- as.vector(unique(ord_tab$sites))
  window_area_plot <- as.data.frame(groupNames)
  
  #  find the convex hulls
  for (i in groupNames) {  #loop among groups
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      convex_hull_pts<-chull(keep) #find row numbers of the hull points
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
      assign(paste0("CHxy_",i),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
    }
  }
  
  all_keep <- ord_tab[,2:3]
  convex_hull_pts <- chull(all_keep)
  convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
  assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
  
  # convert the convex hulls into spatstat "polygonal windows" 
  for (i in groupNames) {
    keep<-ord_tab[ord_tab[,1]==i,2:3] #keep data for one population
    if (length(keep$Axis.1) > 2) {
      CH<-get(paste0("CHxy_",i))
      assign(paste0("Win_",gsub(" ","",i)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      print(paste0(i))
      print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
      
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      
      window_area_plot[window_area_plot$groupNames==i, "window_area"] <- WINDOWAREA
      
      CH<-get(paste0("CHxy_all"))
      assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
      WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
      window_area_plot[window_area_plot$groupNames==i, "total_area"] <- WINDOWAREA
    }
  }
  window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
  write.table(window_area_plot, "ConvexHullAnalysis_sites_PCoA_TAXAbased_relabund_filtsamples_filterIncludeToSPECIESONLY.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  }, error=function(e){})
}

if (groupingFlag == TRUE) {
  for (i in 1:numberofGroups) {
    tryCatch({
    groupnum <- paste0("group",i)
    classgroup <- class(sample_data(phylo11)[[groupnum]])
    if (classgroup != "integer" && classgroup != "numeric") {
      message(paste0("group",i," Convex hulls - phylo11 - NMDS"))
      ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0("group",i), NMDS1, NMDS2)
      ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
      ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
      groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        } 
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_NMDS_TAXAbased_relabund_filtsamples_filterIncludeToSPECIESONLY.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      
      message(paste0("group",i," Convex hulls - phylo11 - PCoA"))
      ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0("group",i), Axis.1, Axis.2)
      ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
      ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
      groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        }
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_PCoA_TAXAbased_relabund_filtsamples_filterIncludeToSPECIESONLY.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
    }
    }, error=function(e){})
  }
}

if (chemDataFlag == TRUE) {
  for (i in chem_headers$V1) {
    tryCatch({
    classgroup <- class(sample_data(phylo11)[[i]])
    if (classgroup != "integer" && classgroup != "numeric") {
      message(paste0("Chem ",i," Convex hulls - phylo11 - NMDS"))
      ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0(i), NMDS1, NMDS2)
      ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
      ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
      groupNames <- as.vector(unique(ord_tab[[i]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        } 
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$NMDS1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_NMDS_TAXAbased_relabund_filtsamples_filterIncludeToSPECIESONLY.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      
      message(paste0("Chem ",i," Convex hulls - phylo11 - PCoA"))
      ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, type="samples")
      ord_tab <- ord_plot$data
      ord_tab <- ord_tab %>% select(paste0(i), Axis.1, Axis.2)
      ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
      ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
      groupNames <- as.vector(unique(ord_tab[[i]]))
      window_area_plot <- as.data.frame(groupNames)
      
      #  find the convex hulls
      for (j in groupNames) {  #loop among groups
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          convex_hull_pts<-chull(keep) #find row numbers of the hull points
          convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
          assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
        }
      }
      
      all_keep <- ord_tab[,2:3]
      convex_hull_pts <- chull(all_keep)
      convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
      assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
      
      # convert the convex hulls into spatstat "polygonal windows" 
      for (j in groupNames) {
        keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
        if (length(keep$Axis.1) > 2) {
          CH<-get(paste0("CHxy_",j))
          assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          print(paste0(j))
          print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
          
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          
          window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
          
          CH<-get(paste0("CHxy_all"))
          assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
          WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
          window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
        }
      }
      window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
      write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_PCoA_TAXAbased_relabund_filtsamples_filterIncludeToSPECIESONLY.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
    }
    }, error=function(e){})
  }
}

setwd(paste0(as.character(args[1]), "/05_Ordination/Taxonomy_merge_based/filterInclude_TOSPECIES_only/relative_abundance", sep = ""))
if (replicateFlag == TRUE) {
  if (groupingFlag == TRUE) {
    for (i in 1:numberofGroups) {
      tryCatch({
      groupnum <- paste0("group",i)
      classgroup <- class(sample_data(phylo15a)[[groupnum]])
      if (classgroup != "integer" && classgroup != "numeric") {
        message(paste0("group",i," Convex hulls - phylo15a - NMDS"))
        ord_plot <- plot_ordination(phylo15a, phylo15a_NMDS.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0("group",i), NMDS1, NMDS2)
        ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
        ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
        groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          } 
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_NMDS_TAXAbased_relabund_replicateGroupedSamples_filterIncludeToSPECIESONLY.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
        
        message(paste0("group",i," Convex hulls - phylo15a - PCoA"))
        ord_plot <- plot_ordination(phylo15a, phylo15a_PCoA.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0("group",i), Axis.1, Axis.2)
        ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
        ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
        groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          }
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_PCoA_TAXAbased_relabund_replicateGroupedSamples_filterIncludeToSPECIESONLY.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      }
      }, error=function(e){})
    }
  }
  
  if (chemDataFlag == TRUE) {
    for (i in chem_headers$V1) {
      tryCatch({
      classgroup <- class(sample_data(phylo15a)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        message(paste0("Chem ",i," Convex hulls - phylo15a - NMDS"))
        ord_plot <- plot_ordination(phylo15a, phylo15a_NMDS.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0(i), NMDS1, NMDS2)
        ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
        ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
        groupNames <- as.vector(unique(ord_tab[[i]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          } 
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_NMDS_TAXAbased_relabund_replicateGroupedSamples_filterIncludeToSPECIESONLY.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
        
        message(paste0("Chem ",i," Convex hulls - phylo15a - PCoA"))
        ord_plot <- plot_ordination(phylo15a, phylo15a_PCoA.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0(i), Axis.1, Axis.2)
        ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
        ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
        groupNames <- as.vector(unique(ord_tab[[i]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          }
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_PCoA_TAXAbased_relabund_replicateGroupedSamples_filterIncludeToSPECIESONLY.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      }
      }, error=function(e){})
    }
  }
}

if (sitelabelFlag == TRUE) {
  if (groupingFlag == TRUE) {
    for (i in 1:numberofGroups) {
      tryCatch({
      groupnum <- paste0("group",i)
      classgroup <- class(sample_data(phylo15b)[[groupnum]])
      if (classgroup != "integer" && classgroup != "numeric") {
        message(paste0("group",i," Convex hulls - phylo15b - NMDS"))
        ord_plot <- plot_ordination(phylo15b, phylo15b_NMDS.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0("group",i), NMDS1, NMDS2)
        ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
        ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
        groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          } 
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_NMDS_TAXAbased_relabund_siteGroupedSamples_filterIncludeToSPECIESONLY.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
        
        message(paste0("group",i," Convex hulls - phylo15b - PCoA"))
        ord_plot <- plot_ordination(phylo15b, phylo15b_PCoA.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0("group",i), Axis.1, Axis.2)
        ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
        ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
        groupNames <- as.vector(unique(ord_tab[[paste0("group",i)]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          }
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_PCoA_TAXAbased_relabund_siteGroupedSamples_filterIncludeToSPECIESONLY.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      }
      }, error=function(e){})
    }
  }
  
  if (chemDataFlag == TRUE) {
    for (i in chem_headers$V1) {
      tryCatch({
      classgroup <- class(sample_data(phylo15b)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        message(paste0("Chem ",i," Convex hulls - phylo15b - NMDS"))
        ord_plot <- plot_ordination(phylo15b, phylo15b_NMDS.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0(i), NMDS1, NMDS2)
        ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
        ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
        groupNames <- as.vector(unique(ord_tab[[i]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          } 
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$NMDS1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_NMDS_TAXAbased_relabund_siteGroupedSamples_filterIncludeToSPECIESONLY.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
        
        message(paste0("Chem ",i," Convex hulls - phylo15b - PCoA"))
        ord_plot <- plot_ordination(phylo15b, phylo15b_PCoA.ord, type="samples")
        ord_tab <- ord_plot$data
        ord_tab <- ord_tab %>% select(paste0(i), Axis.1, Axis.2)
        ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
        ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
        groupNames <- as.vector(unique(ord_tab[[i]]))
        window_area_plot <- as.data.frame(groupNames)
        
        #  find the convex hulls
        for (j in groupNames) {  #loop among groups
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            convex_hull_pts<-chull(keep) #find row numbers of the hull points
            convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1]) #close the polygon
            assign(paste0("CHxy_",j),convexhull.xy(keep)$bdry[[1]]) #save the coordinates of the convex hull
          }
        }
        
        all_keep <- ord_tab[,2:3]
        convex_hull_pts <- chull(all_keep)
        convex_hull_pts<-c(convex_hull_pts,convex_hull_pts[1])
        assign(paste0("CHxy_all"),convexhull.xy(all_keep)$bdry[[1]])
        
        # convert the convex hulls into spatstat "polygonal windows" 
        for (j in groupNames) {
          keep<-ord_tab[ord_tab[,1]==j,2:3] #keep data for one population
          if (length(keep$Axis.1) > 2) {
            CH<-get(paste0("CHxy_",j))
            assign(paste0("Win_",gsub(" ","",j)),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            print(paste0(j))
            print(summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y))))
            
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            
            window_area_plot[window_area_plot$groupNames==j, "window_area"] <- WINDOWAREA
            
            CH<-get(paste0("CHxy_all"))
            assign(paste0("Win_all"),owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))
            WINDOWAREA <- summary(owin(xrange=range(CH$x),yrange=range(CH$y),poly=list(x=CH$x,y=CH$y)))$area
            window_area_plot[window_area_plot$groupNames==j, "total_area"] <- WINDOWAREA
          }
        }
        window_area_plot$percentOfTotal <- 100* window_area_plot$window_area / window_area_plot$total_area
        write.table(window_area_plot, paste0("ConvexHullAnalysis_chem_",i,"_PCoA_TAXAbased_relabund_siteGroupedSamples_filterIncludeToSPECIESONLY.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      }
      }, error=function(e){})
    }
  }
}

#TAXA-based, read counts (phylo10)
setwd(paste0(as.character(args[1]), "/05_Ordination/Taxonomy_merge_based/filterInclude_TOSPECIES_only/read_count", sep = ""))

tryCatch({
ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                            type="taxa", 
                            color = "Phylum",
                            title=paste0("NMDS - stress ",round(phylo10_stress, 3)))
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='NMDS_TAXAbased_taxa_readcount_filtsamples_alltaxa_filterTOSPECIES_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='NMDS_TAXAbased_taxa_readcount_filtsamples_alltaxa_filterTOSPECIES.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                            type="taxa", 
                            color = "Genus",
                            title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
  facet_wrap(~Phylum)
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='NMDS_TAXAbased_taxa_readcount_filtsamples_alltaxa_facetPhylum_colorGenus_filterTOSPECIES_legend.pdf', width = 20, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='NMDS_TAXAbased_taxa_readcount_filtsamples_alltaxa_facetPhylum_colorGenus_filterTOSPECIES.pdf', width = 15, height = 15)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                            type="taxa", 
                            color = "Phylum",
                            title="PCoA")
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='PCoA_TAXAbased_taxa_readcount_filtsamples_alltaxa_filterTOSPECIES_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='PCoA_TAXAbased_taxa_readcount_filtsamples_alltaxa_filterTOSPECIES.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                            type="taxa", 
                            color = "Genus",
                            title="PCoA") +
  facet_wrap(~Phylum)
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='PCoA_TAXAbased_taxa_readcount_filtsamples_alltaxa_facetPhylum_colorGenus_filterTOSPECIES_legend.pdf', width = 20, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='PCoA_TAXAbased_taxa_readcount_filtsamples_alltaxa_facetPhylum_colorGenus_filterTOSPECIES.pdf', width = 15, height = 15)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                            type="samples", 
                            title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
  geom_point(size=3) +
  geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
pdf(file='NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa_filterTOSPECIES.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                            type="samples", 
                            title="PCoA") +
  geom_point(size=3) +
  geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
pdf(file='PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa_filterTOSPECIES.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                            type="biplot", 
                            color = "Phylum",
                            title=paste0("NMDS - stress ",round(phylo10_stress, 3)))
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='NMDS_TAXAbased_taxa_readcount_filtsamples_alltaxa_biplot_filterTOSPECIES_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='NMDS_TAXAbased_taxa_readcount_filtsamples_alltaxa_biplot_filterTOSPECIES.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                            type="biplot", 
                            color = "Phylum",
                            title="PCoA")
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='PCoA_TAXAbased_taxa_readcount_filtsamples_alltaxa_biplot_filterTOSPECIES_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='PCoA_TAXAbased_taxa_readcount_filtsamples_alltaxa_biplot_filterTOSPECIES.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

if (groupingFlag == TRUE) {
  for (i in 1:numberofGroups) {
    tryCatch({
    if (replicateFlag == TRUE) {
      ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored_encircleReplicates_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored_encircleReplicates_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    if (sitelabelFlag == TRUE) {
      ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored_encircleSites_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored_encircleSites_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    groupnum <- paste0("group",i)
    classgroup <- class(sample_data(phylo10)[[groupnum]])
    if (classgroup == "integer" || classgroup == "numeric") {
      ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
        geom_point(size=3)
      pdf(file=paste0('NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3)
      pdf(file=paste0('PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    } else {
      ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group = eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        guides(fill = FALSE)
      pdf(file=paste0('NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored_encircleGroup_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group = eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        guides(fill = FALSE)
      pdf(file=paste0('PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa_group',i,'Colored_encircleGroup_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
  }
}

if (chemDataFlag == TRUE) {
  for (i in chem_headers$V1) {
    tryCatch({
    if (replicateFlag == TRUE) {
      classgroup <- class(sample_data(phylo10)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                                  type="samples",
                                  color=paste0(i),
                                  title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa_chem_',i,'_Colored_encircleReplicates_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                                  type="samples", 
                                  color=paste0(i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa_chem_',i,'_Colored_encircleReplicates_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    if (sitelabelFlag == TRUE) {
      classgroup <- class(sample_data(phylo10)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      
      ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                                  type="samples",
                                  color=paste0(i),
                                  title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa_chem_',i,'_Colored_encircleSites_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                                  type="samples", 
                                  color=paste0(i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa_chem_',i,'_Colored_encircleSites_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    classgroup <- class(sample_data(phylo10)[[i]])
    if (classgroup != "integer" && classgroup != "numeric") {
      viridis_scale <- 'scale_color_viridis_d(option="C")'
    } else {
      viridis_scale <- 'scale_color_viridis_c(option="C")'
    }
    if (classgroup == "integer" || classgroup == "numeric") {
      ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                                  type="samples",
                                  color=paste0(i),
                                  title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
        geom_point(size=3) + eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa_chem_',i,'_Colored_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                                  type="samples", 
                                  color=paste0(i),
                                  title="PCoA") +
        geom_point(size=3) + eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa_chem_',i,'_Colored_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
  }
}

tryCatch({
if (replicateFlag == TRUE) {
  ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                              type="samples",
                              color="replicates",
                              title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa_encircleReplicates_filterTOSPECIES.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  
  ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                              type="samples", 
                              color="replicates",
                              title="PCoA") +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa_encircleReplicates_filterTOSPECIES.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
}
}, error=function(e){})

tryCatch({
if (sitelabelFlag == TRUE) {
  ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                              type="samples",
                              color="sites",
                              title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
    geom_point(size=3) +
    geom_encircle(aes(fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa_encircleSites_filterTOSPECIES.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  
  ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                              type="samples", 
                              color="sites",
                              title="PCoA") +
    geom_point(size=3) +
    geom_encircle(aes(fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa_encircleSites_filterTOSPECIES.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
}
}, error=function(e){})

tryCatch({
if (replicateFlag == TRUE && sitelabelFlag == TRUE) {
  ord_plot <- plot_ordination(phylo10, phylo10_NMDS.ord, 
                              type="samples",
                              color="sites",
                              shape="replicates",
                              title=paste0("NMDS - stress ",round(phylo10_stress, 3))) +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='NMDS_TAXAbased_samples_readcount_filtsamples_alltaxa_encircleshapeReplicates_colorSites_filterTOSPECIES.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  
  ord_plot <- plot_ordination(phylo10, phylo10_PCoA.ord, 
                              type="samples", 
                              color="sites",
                              shape="replicates",
                              title="PCoA") +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='PCoA_TAXAbased_samples_readcount_filtsamples_alltaxa_encircleshapeReplicates_colorSites_filterTOSPECIES.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
}
}, error=function(e){})


if (replicateFlag == TRUE) {
  setwd(paste0(as.character(args[1]), "/05_Ordination/Taxonomy_merge_based/filterInclude_TOSPECIES_only/relative_abundance", sep = ""))
  
  tryCatch({
  ord_plot <- plot_ordination(phylo15a, phylo15a_NMDS.ord, 
                              type="taxa", 
                              color = "Phylum",
                              title=paste0("NMDS - stress ",round(phylo15a_stress, 3)))
  ord_plot_legend <- get_legend(ord_plot)
  ord_plot <- ord_plot + theme(legend.position='none')
  pdf(file='NMDS_TAXAbased_taxa_relabund_replicateGroupedSamples_alltaxa_filterTOSPECIES_legend.pdf', width = 11, height = 6)
  print(as_ggplot(ord_plot_legend))
  dev.off()
  pdf(file='NMDS_TAXAbased_taxa_relabund_replicateGroupedSamples_alltaxa_filterTOSPECIES.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  ord_plot <- plot_ordination(phylo15a, phylo15a_PCoA.ord, 
                              type="taxa", 
                              color = "Phylum",
                              title="PCoA")
  ord_plot_legend <- get_legend(ord_plot)
  ord_plot <- ord_plot + theme(legend.position='none')
  pdf(file='PCoA_TAXAbased_taxa_relabund_replicateGroupedSamples_alltaxa_filterTOSPECIES_legend.pdf', width = 11, height = 6)
  print(as_ggplot(ord_plot_legend))
  dev.off()
  pdf(file='PCoA_TAXAbased_taxa_relabund_replicateGroupedSamples_alltaxa_filterTOSPECIES.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  ord_plot <- plot_ordination(phylo15a, phylo15a_NMDS.ord, 
                              type="samples", 
                              title=paste0("NMDS - stress ",round(phylo15a_stress, 3))) +
    geom_point(size=3) +
    geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
  pdf(file='NMDS_TAXAbased_samples_relabund_replicateGroupedSamples_alltaxa_filterTOSPECIES.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  ord_plot <- plot_ordination(phylo15a, phylo15a_PCoA.ord, 
                              type="samples", 
                              title="PCoA") +
    geom_point(size=3) +
    geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
  pdf(file='PCoA_TAXAbased_samples_relabund_replicateGroupedSamples_alltaxa_filterTOSPECIES.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  if (groupingFlag == TRUE) {
    for (i in 1:numberofGroups) {
      tryCatch({
      groupnum <- paste0("group",i)
      classgroup <- class(sample_data(phylo15a)[[groupnum]])
      if (classgroup == "integer" || classgroup == "numeric") {
        ord_plot <- plot_ordination(phylo15a, phylo15a_NMDS.ord, 
                                    type="samples",
                                    color=paste0("group",i),
                                    title=paste0("NMDS - stress ",round(phylo15a_stress, 3))) +
          geom_point(size=3)
        pdf(file=paste0('NMDS_TAXAbased_samples_relabund_replicateGroupedSamples_alltaxa_group',i,'Colored_filterTOSPECIES.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
        
        ord_plot <- plot_ordination(phylo15a, phylo15a_PCoA.ord, 
                                    type="samples", 
                                    color=paste0("group",i),
                                    title="PCoA") +
          geom_point(size=3)
        pdf(file=paste0('PCoA_TAXAbased_samples_relabund_replicateGroupedSamples_alltaxa_group',i,'Colored_filterTOSPECIES.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
      } else {
        ord_plot <- plot_ordination(phylo15a, phylo15a_NMDS.ord, 
                                    type="samples",
                                    color=paste0("group",i),
                                    title=paste0("NMDS - stress ",round(phylo15a_stress, 3))) +
          geom_point(size=3) +
          geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group = eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
          guides(fill = FALSE)
        pdf(file=paste0('NMDS_TAXAbased_samples_relabund_replicateGroupedSamples_alltaxa_group',i,'Colored_encircleGroup_filterTOSPECIES.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
        
        ord_plot <- plot_ordination(phylo15a, phylo15a_PCoA.ord, 
                                    type="samples", 
                                    color=paste0("group",i),
                                    title="PCoA") +
          geom_point(size=3) +
          geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group = eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
          guides(fill = FALSE)
        pdf(file=paste0('PCoA_TAXAbased_samples_relabund_replicateGroupedSamples_alltaxa_group',i,'Colored_encircleGroup_filterTOSPECIES.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
      }
      }, error=function(e){})
    }
  }
  
  if (chemDataFlag == TRUE) {
    for (i in chem_headers$V1) {
      tryCatch({
      classgroup <- class(sample_data(phylo15a)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      if (classgroup == "integer" || classgroup == "numeric") {
        ord_plot <- plot_ordination(phylo15a, phylo15a_NMDS.ord, 
                                    type="samples",
                                    color=paste0(i),
                                    title=paste0("NMDS - stress ",round(phylo15a_stress, 3))) +
          geom_point(size=3) + eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
        pdf(file=paste0('NMDS_TAXAbased_samples_relabund_replicateGroupedSamples_alltaxa_chem_',i,'_Colored_filterTOSPECIES.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
        
        ord_plot <- plot_ordination(phylo15a, phylo15a_PCoA.ord, 
                                    type="samples", 
                                    color=paste0(i),
                                    title="PCoA") +
          geom_point(size=3) + eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
        pdf(file=paste0('PCoA_TAXAbased_samples_relabund_replicateGroupedSamples_alltaxa_chem_',i,'_Colored_filterTOSPECIES.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
      }
      }, error=function(e){})
    }
  }
}

if (sitelabelFlag == TRUE) {
  setwd(paste0(as.character(args[1]), "/05_Ordination/Taxonomy_merge_based/filterInclude_TOSPECIES_only/relative_abundance", sep = ""))
  
  tryCatch({
  ord_plot <- plot_ordination(phylo15b, phylo15b_NMDS.ord, 
                              type="taxa", 
                              color = "Phylum",
                              title=paste0("NMDS - stress ",round(phylo15b_stress, 3)))
  ord_plot_legend <- get_legend(ord_plot)
  ord_plot <- ord_plot + theme(legend.position='none')
  pdf(file='NMDS_TAXAbased_taxa_relabund_siteGroupedSamples_alltaxa_filterTOSPECIES_legend.pdf', width = 11, height = 6)
  print(as_ggplot(ord_plot_legend))
  dev.off()
  pdf(file='NMDS_TAXAbased_taxa_relabund_siteGroupedSamples_alltaxa_filterTOSPECIES.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  ord_plot <- plot_ordination(phylo15b, phylo15b_PCoA.ord, 
                              type="taxa", 
                              color = "Phylum",
                              title="PCoA")
  ord_plot_legend <- get_legend(ord_plot)
  ord_plot <- ord_plot + theme(legend.position='none')
  pdf(file='PCoA_TAXAbased_taxa_relabund_siteGroupedSamples_alltaxa_filterTOSPECIES_legend.pdf', width = 11, height = 6)
  print(as_ggplot(ord_plot_legend))
  dev.off()
  pdf(file='PCoA_TAXAbased_taxa_relabund_siteGroupedSamples_alltaxa_filterTOSPECIES.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  ord_plot <- plot_ordination(phylo15b, phylo15b_NMDS.ord, 
                              type="samples", 
                              title=paste0("NMDS - stress ",round(phylo15b_stress, 3))) +
    geom_point(size=3) +
    geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
  pdf(file='NMDS_TAXAbased_samples_relabund_siteGroupedSamples_alltaxa_filterTOSPECIES.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  ord_plot <- plot_ordination(phylo15b, phylo15b_PCoA.ord, 
                              type="samples", 
                              title="PCoA") +
    geom_point(size=3) +
    geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
  pdf(file='PCoA_TAXAbased_samples_relabund_siteGroupedSamples_alltaxa_filterTOSPECIES.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  if (groupingFlag == TRUE) {
    for (i in 1:numberofGroups) {
      tryCatch({
      groupnum <- paste0("group",i)
      classgroup <- class(sample_data(phylo15b)[[groupnum]])
      if (classgroup == "integer" || classgroup == "numeric") {
        ord_plot <- plot_ordination(phylo15b, phylo15b_NMDS.ord, 
                                    type="samples",
                                    color=paste0("group",i),
                                    title=paste0("NMDS - stress ",round(phylo15b_stress, 3))) +
          geom_point(size=3)
        pdf(file=paste0('NMDS_TAXAbased_samples_relabund_siteGroupedSamples_alltaxa_group',i,'Colored_filterTOSPECIES.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
        
        ord_plot <- plot_ordination(phylo15b, phylo15b_PCoA.ord, 
                                    type="samples", 
                                    color=paste0("group",i),
                                    title="PCoA") +
          geom_point(size=3)
        pdf(file=paste0('PCoA_TAXAbased_samples_relabund_siteGroupedSamples_alltaxa_group',i,'Colored_filterTOSPECIES.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
      } else {
        ord_plot <- plot_ordination(phylo15b, phylo15b_NMDS.ord, 
                                    type="samples",
                                    color=paste0("group",i),
                                    title=paste0("NMDS - stress ",round(phylo15b_stress, 3))) +
          geom_point(size=3) +
          geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group = eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
          guides(fill = FALSE)
        pdf(file=paste0('NMDS_TAXAbased_samples_relabund_siteGroupedSamples_alltaxa_group',i,'Colored_encircleGroup_filterTOSPECIES.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
        
        ord_plot <- plot_ordination(phylo15b, phylo15b_PCoA.ord, 
                                    type="samples", 
                                    color=paste0("group",i),
                                    title="PCoA") +
          geom_point(size=3) +
          geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group = eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
          guides(fill = FALSE)
        pdf(file=paste0('PCoA_TAXAbased_samples_relabund_siteGroupedSamples_alltaxa_group',i,'Colored_encircleGroup_filterTOSPECIES.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
      }
      }, error=function(e){})
    }
  }
  
  if (chemDataFlag == TRUE) {
    for (i in chem_headers$V1) {
      tryCatch({
      classgroup <- class(sample_data(phylo15b)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      if (classgroup == "integer" || classgroup == "numeric") {
        ord_plot <- plot_ordination(phylo15b, phylo15b_NMDS.ord, 
                                    type="samples",
                                    color=paste0(i),
                                    title=paste0("NMDS - stress ",round(phylo15b_stress, 3))) +
          geom_point(size=3) + eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
        pdf(file=paste0('NMDS_TAXAbased_samples_relabund_siteGroupedSamples_alltaxa_chem_',i,'_Colored_filterTOSPECIES.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
        
        ord_plot <- plot_ordination(phylo15b, phylo15b_PCoA.ord, 
                                    type="samples", 
                                    color=paste0(i),
                                    title="PCoA") +
          geom_point(size=3) + eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
        pdf(file=paste0('PCoA_TAXAbased_samples_relabund_siteGroupedSamples_alltaxa_chem_',i,'_Colored_filterTOSPECIES.pdf'), width = 8, height = 8)
        print(ord_plot)
        dev.off()
      }
      }, error=function(e){})
    }
  }
}

#Taxa-based relative abundance (phylo11)
setwd(paste0(as.character(args[1]), "/05_Ordination/Taxonomy_merge_based/filterInclude_TOSPECIES_only/relative_abundance", sep = ""))

tryCatch({
ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                            type="taxa", 
                            color = "Phylum",
                            title=paste0("NMDS - stress ",round(phylo11_stress, 3)))
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='NMDS_TAXAbased_taxa_relabund_filtsamples_alltaxa_filterTOSPECIES_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='NMDS_TAXAbased_taxa_relabund_filtsamples_alltaxa_filterTOSPECIES.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                            type="taxa", 
                            color = "Genus",
                            title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
  facet_wrap(~Phylum)
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='NMDS_TAXAbased_taxa_relabund_filtsamples_alltaxa_facetPhylum_colorGenus_filterTOSPECIES_legend.pdf', width = 20, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='NMDS_TAXAbased_taxa_relabund_filtsamples_alltaxa_facetPhylum_colorGenus_filterTOSPECIES.pdf', width = 15, height = 15)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                            type="taxa", 
                            color = "Phylum",
                            title="PCoA")
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='PCoA_TAXAbased_taxa_relabund_filtsamples_alltaxa_filterTOSPECIES_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='PCoA_TAXAbased_taxa_relabund_filtsamples_alltaxa_filterTOSPECIES.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                            type="taxa", 
                            color = "Genus",
                            title="PCoA") +
  facet_wrap(~Phylum)
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='PCoA_TAXAbased_taxa_relabund_filtsamples_alltaxa_facetPhylum_colorGenus_filterTOSPECIES_legend.pdf', width = 20, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='PCoA_TAXAbased_taxa_relabund_filtsamples_alltaxa_facetPhylum_colorGenus_filterTOSPECIES.pdf', width = 15, height = 15)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                            type="samples", 
                            title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
  geom_point(size=3) +
  geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
pdf(file='NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa_filterTOSPECIES.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                            type="samples", 
                            title="PCoA") +
  geom_point(size=3) +
  geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
pdf(file='PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa_filterTOSPECIES.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                            type="biplot", 
                            color = "Phylum",
                            title=paste0("NMDS - stress ",round(phylo11_stress, 3)))
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='NMDS_TAXAbased_taxa_relabund_filtsamples_alltaxa_biplot_filterTOSPECIES_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='NMDS_TAXAbased_taxa_relabund_filtsamples_alltaxa_biplot_filterTOSPECIES.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                            type="biplot", 
                            color = "Phylum",
                            title="PCoA")
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='PCoA_TAXAbased_taxa_relabund_filtsamples_alltaxa_biplot_filterTOSPECIES_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='PCoA_TAXAbased_taxa_relabund_filtsamples_alltaxa_biplot_filterTOSPECIES.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

if (groupingFlag == TRUE) {
  for (i in 1:numberofGroups) {
    tryCatch({
    if (replicateFlag == TRUE) {
      ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored_encircleReplicates_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored_encircleReplicates_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    if (sitelabelFlag == TRUE) {
      ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored_encircleSites_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored_encircleSites_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    groupnum <- paste0("group",i)
    classgroup <- class(sample_data(phylo11)[[groupnum]])
    if (classgroup == "integer" || classgroup == "numeric") {
      ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
        geom_point(size=3)
      pdf(file=paste0('NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3)
      pdf(file=paste0('PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    } else {
      ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                                  type="samples",
                                  color=paste0("group",i),
                                  title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group = eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        guides(fill = FALSE)
      pdf(file=paste0('NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored_encircleGroup_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group = eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        guides(fill = FALSE)
      pdf(file=paste0('PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa_group',i,'Colored_encircleGroup_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
  }
}

if (chemDataFlag == TRUE) {
  for (i in chem_headers$V1) {
    tryCatch({
    if (replicateFlag == TRUE) {
      classgroup <- class(sample_data(phylo11)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                                  type="samples",
                                  color=paste0(i),
                                  title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa_chem_',i,'_Colored_encircleReplicates_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                                  type="samples", 
                                  color=paste0(i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa_chem_',i,'_Colored_encircleReplicates_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    if (sitelabelFlag == TRUE) {
      classgroup <- class(sample_data(phylo11)[[i]])
      if (classgroup != "integer" && classgroup != "numeric") {
        viridis_scale <- 'scale_color_viridis_d(option="C")'
      } else {
        viridis_scale <- 'scale_color_viridis_c(option="C")'
      }
      ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                                  type="samples",
                                  color=paste0(i),
                                  title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa_chem_',i,'_Colored_encircleSites_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                                  type="samples", 
                                  color=paste0(i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa_chem_',i,'_Colored_encircleSites_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
    
    tryCatch({
    classgroup <- class(sample_data(phylo11)[[i]])
    if (classgroup != "integer" && classgroup != "numeric") {
      viridis_scale <- 'scale_color_viridis_d(option="C")'
    } else {
      viridis_scale <- 'scale_color_viridis_c(option="C")'
    }
    if (classgroup == "integer" || classgroup == "numeric") {
      ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                                  type="samples",
                                  color=paste0(i),
                                  title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
        geom_point(size=3) + eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa_chem_',i,'_Colored_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                                  type="samples", 
                                  color=paste0(i),
                                  title="PCoA") +
        geom_point(size=3) + eval(parse(text = viridis_scale)) + geom_point(size = 3, shape=1, color = "black")
      pdf(file=paste0('PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa_chem_',i,'_Colored_filterTOSPECIES.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
    }
    }, error=function(e){})
  }
}

tryCatch({
if (replicateFlag == TRUE) {
  ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                              type="samples",
                              color="replicates",
                              title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa_encircleReplicates_filterTOSPECIES.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  
  ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                              type="samples", 
                              color="replicates",
                              title="PCoA") +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa_encircleReplicates_filterTOSPECIES.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
}
}, error=function(e){})

tryCatch({
if (sitelabelFlag == TRUE) {
  ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                              type="samples",
                              color="sites",
                              title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
    geom_point(size=3) +
    geom_encircle(aes(fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa_encircleSites_filterTOSPECIES.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  
  ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                              type="samples", 
                              color="sites",
                              title="PCoA") +
    geom_point(size=3) +
    geom_encircle(aes(fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa_encircleSites_filterTOSPECIES.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
}
}, error=function(e){})
  
tryCatch({
if (replicateFlag == TRUE && sitelabelFlag == TRUE) {
  ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                              type="samples",
                              color="sites",
                              shape="replicates",
                              title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='NMDS_TAXAbased_samples_relabund_filtsamples_alltaxa_encircleshapeReplicates_colorSites_filterTOSPECIES.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  
  ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                              type="samples", 
                              color="sites",
                              shape="replicates",
                              title="PCoA") +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='PCoA_TAXAbased_samples_relabund_filtsamples_alltaxa_encircleshapeReplicates_colorSites_filterTOSPECIES.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
}
}, error=function(e){})




if (original_replicateFlag == TRUE) {
  replicateFlag <- TRUE
}
if (original_sitelabelFlag == TRUE) {
  sitelabelFlag <- TRUE
}

########################
#
# USER SPECIFIED TAXA OF INTEREST
#
########################

##################################
#
#  REMAKE phylo10, phylo11, phylo15a/b
#
##################################
asv_taxonomy <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/ASV2Taxonomy/",as.character(args[3]),"_asvTaxonomyTable_NOUNKNOWNS.txt", sep = ""), header=TRUE, stringsAsFactors=FALSE)
row.names(asv_taxonomy) <- asv_taxonomy$ASV
asv_taxonomy <- asv_taxonomy %>% select(-ASV)
asv_taxonomy_mat <- as.matrix(asv_taxonomy)

sample_metadata <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/sample_metadata_forR.txt", sep = ""), header=TRUE, stringsAsFactors=TRUE)
row.names(sample_metadata) <- sample_metadata$Sample
sample_metadata <- sample_metadata %>% select(-Sample)

TAX = tax_table(asv_taxonomy_mat)
samples = sample_data(sample_metadata)

###################################
##
##  Import Taxonomy-based, raw reads, qual filtered (phylo10)
##
###################################
if (controlFlag == TRUE && filteredLowQualSamples == TRUE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_lowEffortSamplesRemoved.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo10 <- phyloseq(ASV, TAX, samples)
} else if (controlFlag == TRUE && filteredLowQualSamples == FALSE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_controlsRemoved.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo10 <- phyloseq(ASV, TAX, samples)
} else if (controlFlag == FALSE && filteredLowQualSamples == TRUE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_lowEffortSamplesRemoved.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo10 <- phyloseq(ASV, TAX, samples)
} else if (controlFlag == FALSE && filteredLowQualSamples == FALSE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy_NOUNKNOWNS.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo10 <- phyloseq(ASV, TAX, samples)
}

##################################
#
#  Import Taxonomy-based, relative abund, qual filtered (phylo11)
#
##################################
asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
row.names(asv_count) <- asv_count$x
asv_count <- asv_count %>% select(-x)
asv_count_mat <- as.matrix(asv_count)
ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
phylo11 <- phyloseq(ASV, TAX, samples)

##################################
#
#  Import Taxonomy-based, relative abund, qual filtered, grouped on replicates/sites (phylo15a/b)
#
##################################
if (replicateFlag == TRUE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund_groupedByReplicates.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo15a <- phyloseq(ASV, TAX, grouped_samplesA)
}
if (sitelabelFlag == TRUE) {
  asv_count <- read.delim(paste0(as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund_groupedBySites.tsv", sep = ""), header=TRUE, stringsAsFactors=FALSE)
  row.names(asv_count) <- asv_count$x
  asv_count <- asv_count %>% select(-x)
  asv_count_mat <- as.matrix(asv_count)
  ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
  phylo15b <- phyloseq(ASV, TAX, grouped_samplesB)
}



##################################
#
#  SUBSET phylo objects to user interest
#
##################################
if (taxaofinterestFlag == TRUE) {
  taxa_of_interest <- read.delim(as.character(args[14]), header=FALSE, stringsAsFactors=FALSE)
  taxa_of_interest <- as.vector(taxa_of_interest)
  hierarchy_level <- as.character(args[13])
  
  skip_filtered<-TRUE
  
  tryCatch({
    phylo4 <- subset_taxa(phylo4, eval(parse(text = hierarchy_level)) %in% taxa_of_interest$V1)
    phylo5 <- subset_taxa(phylo5, eval(parse(text = hierarchy_level)) %in% taxa_of_interest$V1)
    skip_filtered<-FALSE
  }, error=function(e){})
  
  tryCatch({
  phylo1 <- subset_taxa(phylo1, eval(parse(text = hierarchy_level)) %in% taxa_of_interest$V1)
  }, error=function(e){})
  
  tryCatch({
  phylo2 <- subset_taxa(phylo2, eval(parse(text = hierarchy_level)) %in% taxa_of_interest$V1)
  }, error=function(e){})
  
  tryCatch({
  phylo3 <- subset_taxa(phylo3, eval(parse(text = hierarchy_level)) %in% taxa_of_interest$V1)
  }, error=function(e){})
  
  tryCatch({
  phylo10 <- subset_taxa(phylo10, eval(parse(text = hierarchy_level)) %in% taxa_of_interest$V1)
  }, error=function(e){})
  
  tryCatch({
  phylo11 <- subset_taxa(phylo11, eval(parse(text = hierarchy_level)) %in% taxa_of_interest$V1)
  }, error=function(e){})
  
  tryCatch({
  if (skip_filtered == FALSE) {
  phylo12 <- subset_taxa(phylo12, eval(parse(text = hierarchy_level)) %in% taxa_of_interest$V1)
  phylo13 <- subset_taxa(phylo13, eval(parse(text = hierarchy_level)) %in% taxa_of_interest$V1)
  }
  }, error=function(e){})
  
  
  if (replicateFlag == TRUE) {
    tryCatch({
    phylo7a <- subset_taxa(phylo7a, eval(parse(text = hierarchy_level)) %in% taxa_of_interest$V1)
    }, error=function(e){})
    
    tryCatch({
    if (skip_filtered == FALSE) {
    phylo9a <- subset_taxa(phylo9a, eval(parse(text = hierarchy_level)) %in% taxa_of_interest$V1)
    }
    }, error=function(e){})
    
    tryCatch({
    phylo15a <- subset_taxa(phylo15a, eval(parse(text = hierarchy_level)) %in% taxa_of_interest$V1)
    }, error=function(e){})
    
    tryCatch({
    if (skip_filtered == FALSE) {
    phylo17a <- subset_taxa(phylo17a, eval(parse(text = hierarchy_level)) %in% taxa_of_interest$V1)
    }
    }, error=function(e){})
  }
  
  if (sitelabelFlag == TRUE) {
    tryCatch({
    phylo7b <- subset_taxa(phylo7b, eval(parse(text = hierarchy_level)) %in% taxa_of_interest$V1)
    }, error=function(e){})
    
    tryCatch({
    if (skip_filtered == FALSE) {
    phylo9b <- subset_taxa(phylo9b, eval(parse(text = hierarchy_level)) %in% taxa_of_interest$V1)
    }
    }, error=function(e){})
    
    tryCatch({
    phylo15b <- subset_taxa(phylo15b, eval(parse(text = hierarchy_level)) %in% taxa_of_interest$V1)
    }, error=function(e){})
    
    tryCatch({
    if (skip_filtered == FALSE) {
    phylo17b <- subset_taxa(phylo17b, eval(parse(text = hierarchy_level)) %in% taxa_of_interest$V1)
    }
    }, error=function(e){})
  }
  
  ########################
  #
  # BARCHARTS
  #
  ########################
  setwd(paste0(as.character(args[1]), "/Taxa_of_interest/02_Barcharts/read_count", sep = ""))
  
  tryCatch({
  for (i in c("Phylum")) {
    if (NAremoveFlag == TRUE) {
      bar_plot <- plot_bar(subset_taxa(phylo1, eval(parse(text = paste0("!is.na(",i,")")))), fill = i) +
        geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
        ylab("Abundance (read count)") +
        theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
    } else {
      bar_plot <- plot_bar(phylo1, fill = i) +
        geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
        ylab("Abundance (read count)") +
        theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
    }
    bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = sample_order$V1)
    bar_plot_legend <- get_legend(bar_plot)
    bar_plot <- bar_plot + theme(legend.position='none')
    legend_call <- eval(parse(text = paste0("base_",i,"_count")))
    w_choice <- 6.4 + ((legend_call*2.3)/20)
    pdf(file=paste0('barplot_readcount_allsamples_alltaxa_',i,'_legend.pdf', sep = ""), width = w_choice, height = 6)
    print(as_ggplot(bar_plot_legend))
    dev.off()
    pdf(file=paste0('barplot_readcount_allsamples_alltaxa_',i,'.pdf', sep = ""), width = 11, height = 8.5)
    print(bar_plot)
    dev.off()
  }
  }, error=function(e){})
  
  
  if (controlFlag == TRUE || filteredLowQualSamples == TRUE) {
    tryCatch({
    for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
      if (NAremoveFlag == TRUE) {
        bar_plot <- plot_bar(subset_taxa(phylo2, eval(parse(text = paste0("!is.na(",i,")")))), fill = i) +
          geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
          ylab("Abundance (read count)") +
          theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
      } else {
        bar_plot <- plot_bar(phylo2, fill = i) +
          geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
          ylab("Abundance (read count)") +
          theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
      }
      bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = sample_order$V1)
      bar_plot_legend <- get_legend(bar_plot)
      bar_plot <- bar_plot + theme(legend.position='none')
      legend_call <- eval(parse(text = paste0("base_",i,"_count")))
      w_choice <- 6.4 + ((legend_call*2.3)/20)
      pdf(file=paste0('barplot_readcount_filtsamples_alltaxa_',i,'_legend.pdf', sep = ""), width = w_choice, height = 6)
      print(as_ggplot(bar_plot_legend))
      dev.off()
      pdf(file=paste0('barplot_readcount_filtsamples_alltaxa_',i,'.pdf', sep = ""), width = 11, height = 8.5)
      print(bar_plot)
      dev.off()
    }
    }, error=function(e){})
  } else {
    tryCatch({
    for (i in c("Class", "Order", "Family", "Genus", "Species")) {
      if (NAremoveFlag == TRUE) {
        bar_plot <- plot_bar(subset_taxa(phylo1, eval(parse(text = paste0("!is.na(",i,")")))), fill = i) +
          geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
          ylab("Abundance (read count)") +
          theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
      } else {
        bar_plot <- plot_bar(phylo1, fill = i) +
          geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
          ylab("Abundance (read count)") +
          theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
      }
      bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = sample_order$V1)
      bar_plot_legend <- get_legend(bar_plot)
      bar_plot <- bar_plot + theme(legend.position='none')
      legend_call <- eval(parse(text = paste0("base_",i,"_count")))
      w_choice <- 6.4 + ((legend_call*2.3)/20)
      pdf(file=paste0('barplot_readcount_allsamples_alltaxa_',i,'_legend.pdf', sep = ""), width = w_choice, height = 6)
      print(as_ggplot(bar_plot_legend))
      dev.off()
      pdf(file=paste0('barplot_readcount_allsamples_alltaxa_',i,'.pdf', sep = ""), width = 11, height = 8.5)
      print(bar_plot)
      dev.off()
    }
    }, error=function(e){})
  }
  
  #Filtered taxa barcharts
  tryCatch({
  Phylum_count <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f3 | grep -v ","\"","Phylum","\""," | sort | uniq | wc -l", sep = ""),
                                    intern = TRUE))
  Phylum_NA <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f3 | grep ","\"","^NA$","\""," | sort | uniq | wc -l", sep = ""),
                                 intern = TRUE))
  Class_count <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f4 | grep -v ","\"","Class","\""," | sort | uniq | wc -l", sep = ""),
                                   intern = TRUE))
  Class_NA <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f4 | grep ","\"","^NA$","\""," | sort | uniq | wc -l", sep = ""),
                                intern = TRUE))
  Order_count <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f5 | grep -v ","\"","Order","\""," | sort | uniq | wc -l", sep = ""),
                                   intern = TRUE))
  Order_NA <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f5 | grep ","\"","^NA$","\""," | sort | uniq | wc -l", sep = ""),
                                intern = TRUE))
  Family_count <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f6 | grep -v ","\"","Family","\""," | sort | uniq | wc -l", sep = ""),
                                    intern = TRUE))
  Family_NA <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f6 | grep ","\"","^NA$","\""," | sort | uniq | wc -l", sep = ""),
                                 intern = TRUE))
  Genus_count <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f7 | grep -v ","\"","Genus","\""," | sort | uniq | wc -l", sep = ""),
                                   intern = TRUE))
  Genus_NA <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f7 | grep ","\"","^NA$","\""," | sort | uniq | wc -l", sep = ""),
                                intern = TRUE))
  Species_count <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f8 | grep -v ","\"","Species","\""," | sort | uniq | wc -l", sep = ""),
                                     intern = TRUE))
  Species_NA <- as.numeric(system(paste0("cat ",as.character(args[2]),"/",as.character(args[3]),"/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt | cut -f8 | grep ","\"","^NA$","\""," | sort | uniq | wc -l", sep = ""),
                                  intern = TRUE))
  gg_color <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  if (skip_filtered == FALSE) {
    tryCatch({
  for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    count_call <- eval(parse(text = paste0(i,"_count")))
    NA_call <- eval(parse(text = paste0(i,"_NA")))
    temp_col_pallete <- c(gg_color(count_call-NA_call-1), "grey")
    if (NAremoveFlag == TRUE) {
      bar_plot <- plot_bar(subset_taxa(phylo4, eval(parse(text = paste0("!is.na(",i,")")))), fill = i) +
        geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
        ylab("Abundance (read count)") + scale_fill_manual(values = temp_col_pallete, na.value="lightgrey") + 
        scale_color_manual(values = temp_col_pallete, na.value="lightgrey") +
        theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
    } else {
      bar_plot <- plot_bar(phylo4, fill = i) +
        geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
        ylab("Abundance (read count)") + scale_fill_manual(values = temp_col_pallete, na.value="lightgrey") + 
        scale_color_manual(values = temp_col_pallete, na.value="lightgrey") +
        theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
    }
    bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = sample_order$V1)
    bar_plot_legend <- get_legend(bar_plot)
    bar_plot <- bar_plot + theme(legend.position='none')
    pdf(file=paste0('barplot_readcount_filtsamples_filtLT',filter_percent,'PERCtaxa_',i,'_legend.pdf', sep = ""), width = 11, height = 6)
    print(as_ggplot(bar_plot_legend))
    dev.off()
    pdf(file=paste0('barplot_readcount_filtsamples_filtLT',filter_percent,'PERCtaxa_',i,'.pdf', sep = ""), width = 11, height = 8.5)
    print(bar_plot)
    dev.off()
  }
    }, error=function(e){})
  }
  }, error=function(e){})
  
  
  ###Relative abundance
  setwd(paste0(as.character(args[1]), "/Taxa_of_interest/02_Barcharts/relative_abundance", sep = ""))
  standf = function(x) 100*(x / sum(x))
  
  tryCatch({
  for (i in c("Phylum")) {
    if (NAremoveFlag == TRUE) {
      bar_plot <- plot_bar(transform_sample_counts(subset_taxa(phylo3, eval(parse(text = paste0("!is.na(",i,")")))), standf), fill = i) +
        geom_bar(aes(fill=.data[[i]]), stat="identity", position="stack") +
        ylab("Relative Abundance (%)") +
        theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
    } else {
      bar_plot <- plot_bar(transform_sample_counts(phylo3, standf), fill = i) +
        geom_bar(aes(fill=.data[[i]]), stat="identity", position="stack") +
        ylab("Relative Abundance (%)") +
        theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
    }
    bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = sample_order$V1)
    bar_plot <- bar_plot + theme(legend.position='none')
    pdf(file=paste0('barplot_relabund_filtsamples_alltaxa_',i,'_ASVsShown.pdf', sep = ""), width = 11, height = 8.5)
    print(bar_plot)
    dev.off()
  }
  }, error=function(e){})
  
  tryCatch({
  for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    if (NAremoveFlag == TRUE) {
      bar_plot <- plot_bar(transform_sample_counts(subset_taxa(phylo3, eval(parse(text = paste0("!is.na(",i,")")))), standf), fill = i) +
        geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
        ylab("Relative Abundance (%)") +
        theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
    } else {
      bar_plot <- plot_bar(transform_sample_counts(phylo3, standf), fill = i) +
        geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
        ylab("Relative Abundance (%)") +
        theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
    }
    bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = sample_order$V1)
    bar_plot_legend <- get_legend(bar_plot)
    bar_plot <- bar_plot + theme(legend.position='none')
    legend_call <- eval(parse(text = paste0("base_",i,"_count")))
    w_choice <- 6.4 + ((legend_call*2.3)/20)
    pdf(file=paste0('barplot_relabund_filtsamples_alltaxa_',i,'_legend.pdf', sep = ""), width = w_choice, height = 6)
    print(as_ggplot(bar_plot_legend))
    dev.off()
    pdf(file=paste0('barplot_relabund_filtsamples_alltaxa_',i,'.pdf', sep = ""), width = 11, height = 8.5)
    print(bar_plot)
    dev.off()
  }
  }, error=function(e){})
  
  if (skip_filtered == FALSE) {
    tryCatch({
  for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    count_call <- eval(parse(text = paste0(i,"_count")))
    NA_call <- eval(parse(text = paste0(i,"_NA")))
    temp_col_pallete <- c(gg_color(count_call-NA_call-1), "grey")
    if (NAremoveFlag == TRUE) {
      bar_plot <- plot_bar(transform_sample_counts(subset_taxa(phylo5, eval(parse(text = paste0("!is.na(",i,")")))), standf), fill = i) +
        geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
        ylab("Relative Abundance (%)") + scale_fill_manual(values = temp_col_pallete, na.value="lightgrey") + 
        scale_color_manual(values = temp_col_pallete, na.value="lightgrey") +
        theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
    } else {
      bar_plot <- plot_bar(transform_sample_counts(phylo5, standf), fill = i) +
        geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
        ylab("Relative Abundance (%)") + scale_fill_manual(values = temp_col_pallete, na.value="lightgrey") + 
        scale_color_manual(values = temp_col_pallete, na.value="lightgrey") +
        theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
    }
    bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = sample_order$V1)
    bar_plot_legend <- get_legend(bar_plot)
    bar_plot <- bar_plot + theme(legend.position='none')
    pdf(file=paste0('barplot_relabund_filtsamples_filtLT',filter_percent,'PERCtaxa_',i,'_legend.pdf', sep = ""), width = 11, height = 6)
    print(as_ggplot(bar_plot_legend))
    dev.off()
    pdf(file=paste0('barplot_relabund_filtsamples_filtLT',filter_percent,'PERCtaxa_',i,'.pdf', sep = ""), width = 11, height = 8.5)
    print(bar_plot)
    dev.off()
  }
    }, error=function(e){})
  }
  
  #GROUPED
  tryCatch({
  if (replicateFlag == TRUE) {
    for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
      if (NAremoveFlag == TRUE) {
        bar_plot <- plot_bar(transform_sample_counts(subset_taxa(phylo7a, eval(parse(text = paste0("!is.na(",i,")")))), standf), fill = i) +
          geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
          ylab("Relative Abundance (%)") +
          theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
      } else {
        bar_plot <- plot_bar(transform_sample_counts(phylo7a, standf), fill = i) +
          geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
          ylab("Relative Abundance (%)") +
          theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
      }
      bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = replicateorder$`sample_metadata$replicates`)
      bar_plot_legend <- get_legend(bar_plot)
      bar_plot <- bar_plot + theme(legend.position='none')
      legend_call <- eval(parse(text = paste0("base_",i,"_count")))
      w_choice <- 6.4 + ((legend_call*2.3)/20)
      pdf(file=paste0('barplot_relabund_replicateGroupedSamples_alltaxa_',i,'_legend.pdf', sep = ""), width = w_choice, height = 6)
      print(as_ggplot(bar_plot_legend))
      dev.off()
      pdf(file=paste0('barplot_relabund_replicateGroupedSamples_alltaxa_',i,'.pdf', sep = ""), width = 11, height = 8.5)
      print(bar_plot)
      dev.off()
    }
  }
  }, error=function(e){})
  
  tryCatch({
  if (sitelabelFlag == TRUE) {
    for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
      if (NAremoveFlag == TRUE) {
        bar_plot <- plot_bar(transform_sample_counts(subset_taxa(phylo7b, eval(parse(text = paste0("!is.na(",i,")")))), standf), fill = i) +
          geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
          ylab("Relative Abundance (%)") +
          theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
      } else {
        bar_plot <- plot_bar(transform_sample_counts(phylo7b, standf), fill = i) +
          geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
          ylab("Relative Abundance (%)") +
          theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
      }
      bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = siteorder$`sample_metadata$sites`)
      bar_plot_legend <- get_legend(bar_plot)
      bar_plot <- bar_plot + theme(legend.position='none')
      legend_call <- eval(parse(text = paste0("base_",i,"_count")))
      w_choice <- 6.4 + ((legend_call*2.3)/20)
      pdf(file=paste0('barplot_relabund_siteGroupedSamples_alltaxa_',i,'_legend.pdf', sep = ""), width = w_choice, height = 6)
      print(as_ggplot(bar_plot_legend))
      dev.off()
      pdf(file=paste0('barplot_relabund_siteGroupedSamples_alltaxa_',i,'.pdf', sep = ""), width = 11, height = 8.5)
      print(bar_plot)
      dev.off()
    }
  }
  }, error=function(e){})
  
  tryCatch({
  if (replicateFlag == TRUE) {
    if (skip_filtered == FALSE) {
    for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
      count_call <- eval(parse(text = paste0(i,"_count")))
      NA_call <- eval(parse(text = paste0(i,"_NA")))
      temp_col_pallete <- c(gg_color(count_call-NA_call-1), "grey")
      if (NAremoveFlag == TRUE) {
        bar_plot <- plot_bar(transform_sample_counts(subset_taxa(phylo9a, eval(parse(text = paste0("!is.na(",i,")")))), standf), fill = i) +
          geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
          ylab("Relative Abundance (%)") + scale_fill_manual(values = temp_col_pallete, na.value="lightgrey") + 
          scale_color_manual(values = temp_col_pallete, na.value="lightgrey") +
          theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
      } else {
        bar_plot <- plot_bar(transform_sample_counts(phylo9a, standf), fill = i) +
          geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
          ylab("Relative Abundance (%)") + scale_fill_manual(values = temp_col_pallete, na.value="lightgrey") + 
          scale_color_manual(values = temp_col_pallete, na.value="lightgrey") +
          theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
      }
      bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = replicateorder$`sample_metadata$replicates`)
      bar_plot_legend <- get_legend(bar_plot)
      bar_plot <- bar_plot + theme(legend.position='none')
      pdf(file=paste0('barplot_relabund_replicateGroupedSamples_filtLT',filter_percent,'PERCtaxa_',i,'_legend.pdf', sep = ""), width = 11, height = 6)
      print(as_ggplot(bar_plot_legend))
      dev.off()
      pdf(file=paste0('barplot_relabund_replicateGroupedSamples_filtLT',filter_percent,'PERCtaxa_',i,'.pdf', sep = ""), width = 11, height = 8.5)
      print(bar_plot)
      dev.off()
    }
    }
  }
  }, error=function(e){})
  
  tryCatch({
  if (sitelabelFlag == TRUE) {
    if (skip_filtered == FALSE) {
    for (i in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
      count_call <- eval(parse(text = paste0(i,"_count")))
      NA_call <- eval(parse(text = paste0(i,"_NA")))
      temp_col_pallete <- c(gg_color(count_call-NA_call-1), "grey")
      if (NAremoveFlag == TRUE) {
        bar_plot <- plot_bar(transform_sample_counts(subset_taxa(phylo9b, eval(parse(text = paste0("!is.na(",i,")")))), standf), fill = i) +
          geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
          ylab("Relative Abundance (%)") + scale_fill_manual(values = temp_col_pallete, na.value="lightgrey") + 
          scale_color_manual(values = temp_col_pallete, na.value="lightgrey") +
          theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
      } else {
        bar_plot <- plot_bar(transform_sample_counts(phylo9b, standf), fill = i) +
          geom_bar(aes(color=.data[[i]], fill=.data[[i]]), stat="identity", position="stack") +
          ylab("Relative Abundance (%)") + scale_fill_manual(values = temp_col_pallete, na.value="lightgrey") + 
          scale_color_manual(values = temp_col_pallete, na.value="lightgrey") +
          theme(panel.grid.major.x = element_blank(), axis.text.x=element_text(vjust=0.5))
      }
      bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = siteorder$`sample_metadata$sites`)
      bar_plot_legend <- get_legend(bar_plot)
      bar_plot <- bar_plot + theme(legend.position='none')
      pdf(file=paste0('barplot_relabund_siteGroupedSamples_filtLT',filter_percent,'PERCtaxa_',i,'_legend.pdf', sep = ""), width = 11, height = 6)
      print(as_ggplot(bar_plot_legend))
      dev.off()
      pdf(file=paste0('barplot_relabund_siteGroupedSamples_filtLT',filter_percent,'PERCtaxa_',i,'.pdf', sep = ""), width = 11, height = 8.5)
      print(bar_plot)
      dev.off()
    }
    }
  }
  }, error=function(e){})
  
  ########################
  #
  # HEATMAPS
  #
  ########################
  setwd(paste0(as.character(args[1]), "/Taxa_of_interest/03_Heatmaps/ASV_based", sep = ""))
  standf = function(x) 100*(x / sum(x))
  
  tryCatch({
  heatmap <- plot_heatmap(phylo3, method = "NMDS", distance = "jaccard", 
                          low = "ghostwhite", 
                          high = "red", 
                          na.value = "white")
  pdf(file='heatmap_ASVs_relabund_filtsamples_alltaxa_clustSamples.pdf', width = 11, height = 8.5)
  print(heatmap)
  dev.off()
  heatmap$data$Sample <- factor(heatmap$data$Sample, levels = sample_order$V1)
  pdf(file='heatmap_ASVs_relabund_filtsamples_alltaxa_orderedSamples.pdf', width = 11, height = 8.5)
  print(heatmap)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  if (replicateFlag == TRUE) {
    heatmap <- plot_heatmap(phylo7a, method = "NMDS", distance = "jaccard", 
                            low = "ghostwhite", 
                            high = "red", 
                            na.value = "white")
    pdf(file='heatmap_ASVs_relabund_replicateGroupedSamples_alltaxa_clustReplicates.pdf', width = 11, height = 8.5)
    print(heatmap)
    dev.off()
    heatmap$data$Sample <- factor(heatmap$data$Sample, levels = replicateorder$`sample_metadata$replicates`)
    pdf(file='heatmap_ASVs_relabund_replicateGroupedSamples_alltaxa_orderedReplicates.pdf', width = 11, height = 8.5)
    print(heatmap)
    dev.off()
  }
  }, error=function(e){})
  
  tryCatch({
  if (sitelabelFlag == TRUE) {
    heatmap <- plot_heatmap(phylo7b, method = "NMDS", distance = "jaccard", 
                            low = "ghostwhite", 
                            high = "red", 
                            na.value = "white")
    pdf(file='heatmap_ASVs_relabund_siteGroupedSamples_alltaxa_clustSites.pdf', width = 11, height = 8.5)
    print(heatmap)
    dev.off()
    heatmap$data$Sample <- factor(heatmap$data$Sample, levels = siteorder$`sample_metadata$sites`)
    pdf(file='heatmap_ASVs_relabund_siteGroupedSamples_alltaxa_orderedSites.pdf', width = 11, height = 8.5)
    print(heatmap)
    dev.off()
  }
  }, error=function(e){})
  
  
  setwd(paste0(as.character(args[1]), "/Taxa_of_interest/03_Heatmaps/Taxonomy_merge_based", sep = ""))
  
  tryCatch({
  heatmap <- plot_heatmap(phylo11, method = "NMDS", distance = "jaccard", 
                          low = "ghostwhite", 
                          high = "red", 
                          na.value = "white")
  pdf(file='heatmap_Taxa_relabund_filtsamples_alltaxa_clustSamples.pdf', width = 11, height = 8.5)
  print(heatmap)
  dev.off()
  heatmap$data$Sample <- factor(heatmap$data$Sample, levels = sample_order$V1)
  pdf(file='heatmap_Taxa_relabund_filtsamples_alltaxa_orderedSamples.pdf', width = 11, height = 8.5)
  print(heatmap)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  if (replicateFlag == TRUE) {
    heatmap <- plot_heatmap(phylo15a, method = "NMDS", distance = "jaccard", 
                            low = "ghostwhite", 
                            high = "red", 
                            na.value = "white")
    pdf(file='heatmap_Taxa_relabund_replicateGroupedSamples_alltaxa_clustReplicates.pdf', width = 11, height = 8.5)
    print(heatmap)
    dev.off()
    heatmap$data$Sample <- factor(heatmap$data$Sample, levels = replicateorder$`sample_metadata$replicates`)
    pdf(file='heatmap_Taxa_relabund_replicateGroupedSamples_alltaxa_orderedReplicates.pdf', width = 11, height = 8.5)
    print(heatmap)
    dev.off()
  }
  }, error=function(e){})
  
  tryCatch({
  if (sitelabelFlag == TRUE) {
    heatmap <- plot_heatmap(phylo15b, method = "NMDS", distance = "jaccard", 
                            low = "ghostwhite", 
                            high = "red", 
                            na.value = "white")
    pdf(file='heatmap_Taxa_relabund_siteGroupedSamples_alltaxa_clustSites.pdf', width = 11, height = 8.5)
    print(heatmap)
    dev.off()
    heatmap$data$Sample <- factor(heatmap$data$Sample, levels = siteorder$`sample_metadata$sites`)
    pdf(file='heatmap_Taxa_relabund_siteGroupedSamples_alltaxa_orderedSites.pdf', width = 11, height = 8.5)
    print(heatmap)
    dev.off()
  }
  }, error=function(e){})
  
  ########################
  #
  # NETWORK ANALYSIS
  #
  ########################
  setwd(paste0(as.character(args[1]), "/Taxa_of_interest/06_Network/ASV_based/read_count", sep = ""))
  
  for (i in 1:9) {
    tryCatch({
      taxa_net <- make_network(phylo2, type="taxa", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
      net_plot <- plot_network(taxa_net, phylo2)
      pdf(file=paste0('Network_ASVbased_taxa_readcount_filtsamples_alltaxa_dist0.',i,'_ASVLabeled.pdf'), width = 8, height = 8)
      print(net_plot)
      dev.off()
    }, error=function(e){})
  }
  
  for (i in 1:9) {
    tryCatch({
      sample_net <- make_network(phylo2, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
      net_plot <- plot_network(sample_net, phylo2)
      pdf(file=paste0('Network_ASVbased_samples_readcount_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled.pdf'), width = 8, height = 8)
      print(net_plot)
      dev.off()
    }, error=function(e){})
  }
  
  if (replicateFlag == TRUE) {
    for (i in 1:9) {
      tryCatch({
        sample_net <- make_network(phylo2, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
        net_plot <- plot_network(sample_net, phylo2, color="replicates")
        pdf(file=paste0('Network_ASVbased_samples_readcount_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled_colorReplicates.pdf'), width = 8, height = 8)
        print(net_plot)
        dev.off()
      }, error=function(e){})
    }
  }
  
  if (sitelabelFlag == TRUE) {
    for (i in 1:9) {
      tryCatch({
        sample_net <- make_network(phylo2, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
        net_plot <- plot_network(sample_net, phylo2, color="sites")
        pdf(file=paste0('Network_ASVbased_samples_readcount_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled_colorSites.pdf'), width = 8, height = 8)
        print(net_plot)
        dev.off()
      }, error=function(e){})
    }
  }
  
  setwd(paste0(as.character(args[1]), "/Taxa_of_interest/06_Network/ASV_based/relative_abundance", sep = ""))
  
  for (i in 1:9) {
    tryCatch({
      taxa_net <- make_network(filter_taxa(phylo3, function(x) sum(x >= 1) > 0, TRUE), type="taxa", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
      net_plot <- plot_network(taxa_net, filter_taxa(phylo3, function(x) sum(x >= 1) > 0, TRUE))
      pdf(file=paste0('Network_ASVbased_taxa_relabund_filtsamples_gte1perctaxa_dist0.',i,'_ASVLabeled.pdf'), width = 8, height = 8)
      print(net_plot)
      dev.off()
    }, error=function(e){})
  }
  
  for (i in 1:9) {
    tryCatch({
      sample_net <- make_network(phylo3, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
      net_plot <- plot_network(sample_net, phylo3)
      pdf(file=paste0('Network_ASVbased_samples_relabund_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled.pdf'), width = 8, height = 8)
      print(net_plot)
      dev.off()
    }, error=function(e){})
  }
  
  if (replicateFlag == TRUE) {
    for (i in 1:9) {
      tryCatch({
        sample_net <- make_network(phylo3, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
        net_plot <- plot_network(sample_net, phylo3, color="replicates")
        pdf(file=paste0('Network_ASVbased_samples_relabund_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled_colorReplicates.pdf'), width = 8, height = 8)
        print(net_plot)
        dev.off()
      }, error=function(e){})
    }
  }
  
  if (sitelabelFlag == TRUE) {
    for (i in 1:9) {
      tryCatch({
        sample_net <- make_network(phylo3, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
        net_plot <- plot_network(sample_net, phylo3, color="sites")
        pdf(file=paste0('Network_ASVbased_samples_relabund_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled_colorSites.pdf'), width = 8, height = 8)
        print(net_plot)
        dev.off()
      }, error=function(e){})
    }
  }
  
  setwd(paste0(as.character(args[1]), "/Taxa_of_interest/06_Network/Taxonomy_merge_based/read_count", sep = ""))
  
  for (i in 1:9) {
    tryCatch({
      taxa_net <- make_network(filter_taxa(phylo10, function(x) sum(x >= 1) > 0, TRUE), type="taxa", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
      net_plot <- plot_network(taxa_net, filter_taxa(phylo10, function(x) sum(x >= 1) > 0, TRUE))
      pdf(file=paste0('Network_TAXAbased_taxa_readcount_filtsamples_gte1perctaxa_dist0.',i,'_ASVLabeled.pdf'), width = 8, height = 8)
      print(net_plot)
      dev.off()
    }, error=function(e){})
  }
  
  for (i in 1:9) {
    tryCatch({
      sample_net <- make_network(phylo10, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
      net_plot <- plot_network(sample_net, phylo10)
      pdf(file=paste0('Network_TAXAbased_samples_readcount_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled.pdf'), width = 8, height = 8)
      print(net_plot)
      dev.off()
    }, error=function(e){})
  }
  
  if (replicateFlag == TRUE) {
    for (i in 1:9) {
      tryCatch({
        sample_net <- make_network(phylo10, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
        net_plot <- plot_network(sample_net, phylo10, color="replicates")
        pdf(file=paste0('Network_TAXAbased_samples_readcount_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled_colorReplicates.pdf'), width = 8, height = 8)
        print(net_plot)
        dev.off()
      }, error=function(e){})
    }
  }
  
  if (sitelabelFlag == TRUE) {
    for (i in 1:9) {
      tryCatch({
        sample_net <- make_network(phylo10, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
        net_plot <- plot_network(sample_net, phylo10, color="sites")
        pdf(file=paste0('Network_TAXAbased_samples_readcount_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled_colorSites.pdf'), width = 8, height = 8)
        print(net_plot)
        dev.off()
      }, error=function(e){})
    }
  }
  
  setwd(paste0(as.character(args[1]), "/Taxa_of_interest/06_Network/Taxonomy_merge_based/relative_abundance", sep = ""))
  
  for (i in 1:9) {
    tryCatch({
      taxa_net <- make_network(filter_taxa(phylo11, function(x) sum(x >= 1) > 0, TRUE), type="taxa", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
      net_plot <- plot_network(taxa_net, filter_taxa(phylo11, function(x) sum(x >= 1) > 0, TRUE))
      pdf(file=paste0('Network_TAXAbased_taxa_relabund_filtsamples_gte1perctaxa_dist0.',i,'_ASVLabeled.pdf'), width = 8, height = 8)
      print(net_plot)
      dev.off()
    }, error=function(e){})
  }
  
  for (i in 1:9) {
    tryCatch({
      sample_net <- make_network(phylo11, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
      net_plot <- plot_network(sample_net, phylo11)
      pdf(file=paste0('Network_TAXAbased_samples_relabund_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled.pdf'), width = 8, height = 8)
      print(net_plot)
      dev.off()
    }, error=function(e){})
  }
  
  if (replicateFlag == TRUE) {
    for (i in 1:9) {
      tryCatch({
        sample_net <- make_network(phylo11, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
        net_plot <- plot_network(sample_net, phylo11, color="replicates")
        pdf(file=paste0('Network_TAXAbased_samples_relabund_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled_colorReplicates.pdf'), width = 8, height = 8)
        print(net_plot)
        dev.off()
      }, error=function(e){})
    }
  }
  
  if (sitelabelFlag == TRUE) {
    for (i in 1:9) {
      tryCatch({
        sample_net <- make_network(phylo11, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
        net_plot <- plot_network(sample_net, phylo11, color="sites")
        pdf(file=paste0('Network_TAXAbased_samples_relabund_filtsamples_gte1perctaxa_dist0.',i,'_SampleLabeled_colorSites.pdf'), width = 8, height = 8)
        print(net_plot)
        dev.off()
      }, error=function(e){})
    }
  }
}


#OLD METHOD TO FIND HULL
# df <- as.data.frame(phylo_combined.ord$vectors)
# sample_metadata <- tibble::rownames_to_column(sample_metadata)
# df <- tibble::rownames_to_column(df)
# df <- left_join(select(sample_metadata, rowname, replicates), df, by = "rowname" )
# find_hull <- function(df) df[chull(df$Axis.1, df$Axis.2), ]
# hulls <- ddply(df, "replicates", find_hull)
#in ggplot:
# + geom_polygon(data = hulls, aes(x = Axis.1, y = Axis.2, fill=replicates), alpha = 0.5)


