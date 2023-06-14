#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

########################################
# TEMP WHILE WORKING ON SCRIPT
# args[1]<-"/Users/mcallister/Desktop/test_figs" #FIGURE OUT directory
# args[2]<-"/Users/mcallister/Desktop/chris_test/18S/CP_all_out/ASV2Taxonomy/ASVs_counts_NOUNKNOWNS.tsv" #ASV count table
# args[3]<-"/Users/mcallister/Desktop/chris_test/18S/CP_all_out/ASV2Taxonomy/CP_all_out_asvTaxonomyTable_NOUNKNOWNS.txt" #ASV taxonomy table
# args[4]<-"/Users/mcallister/Desktop/chris_test/18S/CP_all_out/sample_metadata_forR.txt" #sample metadata file
# args[5]<-5 #percent filter (ASVs)
# args[6]<-TRUE #controllabelFlag
# args[7]<-TRUE #low quality sample filter
# args[8]<-30 #percent filter (Samples), based on median (excluding controls)
# args[9]<-"/Users/mcallister/Desktop/chris_test/18S/CP_all_out/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy_NOUNKNOWNS.tsv" #ASV counts merged on taxonomy
# args[10]<- FALSE #Sites metadata category
# args[11]<- TRUE #Replicates metadata category
# args[12]<-"/Users/mcallister/Desktop/chris_test/18S/CP_all_out/dada2/ASVs_counts.tsv"
########################################
library("dplyr")

setwd(paste0(as.character(args[1]), "/processed_tables", sep = ""))

asv_count <- read.delim(as.character(args[2]), header=TRUE, stringsAsFactors=FALSE)
asv_count_WITHUNKNOWNS <- read.delim(as.character(args[12]), header=TRUE, stringsAsFactors=FALSE)
asv_taxonomy <- read.delim(as.character(args[3]), header=TRUE, stringsAsFactors=FALSE)
sample_metadata <- read.delim(as.character(args[4]), header=TRUE, stringsAsFactors=FALSE)

#Control filtering
controlFlag <- as.logical(args[6])

if (controlFlag == TRUE) {
  controls_df <- sample_metadata %>% select(Sample, controls)
  controls_df <- controls_df %>% filter(!is.na(controls))
  control_samples <- as.character(c(controls_df[1:nrow(controls_df), 1]))
  
  #TODO: CONTROL FILTERING (TBD)
  
  asv_count <- asv_count %>% select(-all_of(control_samples))
  asv_count$readsum <- rowSums(asv_count[,-1])
  asv_count <- asv_count %>% filter(readsum > 0)
  asv_count <- asv_count %>% select(-readsum)
  asv_count_WITHUNKNOWNS <- asv_count_WITHUNKNOWNS %>% select(-all_of(control_samples))
  asv_count_WITHUNKNOWNS$readsum <- rowSums(asv_count_WITHUNKNOWNS[,-1])
  asv_count_WITHUNKNOWNS <- asv_count_WITHUNKNOWNS %>% filter(readsum > 0)
  asv_count_WITHUNKNOWNS <- asv_count_WITHUNKNOWNS %>% select(-readsum)
  
  write.table(asv_count, "ASVs_counts_NOUNKNOWNS_controlsRemoved.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  write.table(asv_count_WITHUNKNOWNS, "ASVs_counts_controlsRemoved.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}

#Filter underperforming samples
filterLowQualSamplesFlag <- as.logical(args[7])
sampleFilterPerc <- as.numeric(args[8])

asv_count_samp <- asv_count
row.names(asv_count_samp) <- asv_count_samp$x
asv_count_samp <- asv_count_samp %>% select(-x)
asv_count_samp <- as.data.frame(t(asv_count_samp))
asv_count_samp <- tibble::rownames_to_column(asv_count_samp, "x")
asv_count_samp$readSum <- rowSums(asv_count_samp[,-1])
sample_sums <- as.numeric(asv_count_samp$readSum)
median_sample_effort <- median(sample_sums)

asv_count_samp_WU <- asv_count_WITHUNKNOWNS
row.names(asv_count_samp_WU) <- asv_count_samp_WU$x
asv_count_samp_WU <- asv_count_samp_WU %>% select(-x)
asv_count_samp_WU <- as.data.frame(t(asv_count_samp_WU))
asv_count_samp_WU <- tibble::rownames_to_column(asv_count_samp_WU, "x")
asv_count_samp_WU$readSum <- rowSums(asv_count_samp_WU[,-1])
sample_sums_WU <- as.numeric(asv_count_samp_WU$readSum)
median_sample_effort_WU <- median(sample_sums_WU)

if (filterLowQualSamplesFlag == TRUE) {
  asv_count_samp <- asv_count_samp %>% filter(readSum >= median_sample_effort * (sampleFilterPerc/100))
  asv_count_samp <- asv_count_samp %>% select(-readSum)
  row.names(asv_count_samp) <- asv_count_samp$x
  asv_count_samp <- asv_count_samp %>% select(-x)
  asv_count_samp <- as.data.frame(t(asv_count_samp))
  asv_count_samp <- tibble::rownames_to_column(asv_count_samp, "x")
  asv_count <- asv_count_samp
  write.table(asv_count, "ASVs_counts_NOUNKNOWNS_lowEffortSamplesRemoved.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)

  asv_count_samp_WU <- asv_count_samp_WU %>% filter(readSum >= median_sample_effort_WU * (sampleFilterPerc/100))
  asv_count_samp_WU <- asv_count_samp_WU %>% select(-readSum)
  row.names(asv_count_samp_WU) <- asv_count_samp_WU$x
  asv_count_samp_WU <- asv_count_samp_WU %>% select(-x)
  asv_count_samp_WU <- as.data.frame(t(asv_count_samp_WU))
  asv_count_samp_WU <- tibble::rownames_to_column(asv_count_samp_WU, "x")
  asv_count_WITHUNKNOWNS <- asv_count_samp_WU #NOT CURRENTLY PLANNING TO CONVERT THIS TO REL ABUND OR CONTINUE FURTHER. Consider later.
  write.table(asv_count_WITHUNKNOWNS, "ASVs_counts_lowEffortSamplesRemoved.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  
}

#Convert reads to rel. abund. (%)
asv_countrel <- asv_count
row.names(asv_countrel) <- asv_countrel$x
asv_countrel <- asv_countrel %>% select(-x)
asv_countrel <- as.data.frame(t(asv_countrel))
asv_countrel <- tibble::rownames_to_column(asv_countrel, "x")
asv_countrel$readSum <- rowSums(asv_countrel[,-1])
asv_countrel <- asv_countrel %>% mutate_at(vars(-x, -readSum), list(~(100*./readSum)))
asv_countrel <- asv_countrel %>% select(-readSum)
asv_countrel_rotated1 <- asv_countrel
asv_countrel_rotated2 <- asv_countrel
row.names(asv_countrel) <- asv_countrel$x
asv_countrel <- asv_countrel %>% select(-x)
asv_countrel <- as.data.frame(t(asv_countrel))
asv_countrel <- tibble::rownames_to_column(asv_countrel, "x")

write.table(asv_countrel, "ASVs_counts_NOUNKNOWNS_percentabund.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)


#Grouped counts by mean (not sum) of rel. abund
replicateFlag <- as.logical(args[11])
sitelabelFlag <- as.logical(args[10])

if (sitelabelFlag == TRUE) { #Will expect "sites" column
  asv_countrel_rotated <- asv_countrel_rotated1
  asv_countrel_rotated <- asv_countrel_rotated %>% rename(Sample = x)
  asv_countrel_join <- left_join(select(sample_metadata, Sample, sites), asv_countrel_rotated, by = "Sample")
  temp_colnames <- colnames(asv_countrel_join)
  asv_countrel_join <- asv_countrel_join %>% filter(!is.na(eval(parse(text = temp_colnames[3]))))
  asv_countrel_join <- asv_countrel_join %>% select(-Sample)
  rep_group <- asv_countrel_join %>% group_by(sites) %>% 
    summarise_all(mean)
  rep_group <- rep_group %>% rename(x = sites)
  rep_group <- as.data.frame(rep_group)
  row.names(rep_group) <- rep_group$x
  rep_group <- rep_group %>% select(-x)
  rep_group <- as.data.frame(t(rep_group))
  rep_group <- tibble::rownames_to_column(rep_group, "x")
  write.table(rep_group, "ASVs_counts_NOUNKNOWNS_percentabund_groupedBySites.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  
  sample_repgroup <- sample_metadata
  sample_repgroup <- sample_repgroup %>% select(-Sample)
  sample_repgroup <- sample_repgroup %>% filter(!is.na(sites))
  sample_repgroup <- sample_repgroup %>% group_by(sites) %>%
    summarise_all(list(~toString(unique(.)))) %>% rename(Sample = sites)
  write.table(sample_repgroup, "sample_metadata_NOUNKNOWNS_percentabund_groupedBySites.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}

if (replicateFlag == TRUE) { #Will expect "replicates" column
  asv_countrel_rotated <- asv_countrel_rotated2
  asv_countrel_rotated <- asv_countrel_rotated %>% rename(Sample = x)
  asv_countrel_join <- left_join(select(sample_metadata, Sample, replicates), asv_countrel_rotated, by = "Sample")
  temp_colnames <- colnames(asv_countrel_join)
  asv_countrel_join <- asv_countrel_join %>% filter(!is.na(eval(parse(text = temp_colnames[3]))))
  asv_countrel_join <- asv_countrel_join %>% select(-Sample)
  rep_group <- asv_countrel_join %>% group_by(replicates) %>% 
    summarise_all(mean)
  rep_group <- rep_group %>% rename(x = replicates)
  rep_group <- as.data.frame(rep_group)
  row.names(rep_group) <- rep_group$x
  rep_group <- rep_group %>% select(-x)
  rep_group <- as.data.frame(t(rep_group))
  rep_group <- tibble::rownames_to_column(rep_group, "x")
  write.table(rep_group, "ASVs_counts_NOUNKNOWNS_percentabund_groupedByReplicates.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  
  sample_repgroup <- sample_metadata
  sample_repgroup <- sample_repgroup %>% select(-Sample)
  sample_repgroup <- sample_repgroup %>% filter(!is.na(replicates))
  sample_repgroup <- sample_repgroup %>% group_by(replicates) %>%
    summarise_all(list(~toString(unique(.)))) %>% rename(Sample = replicates)
  write.table(sample_repgroup, "sample_metadata_NOUNKNOWNS_percentabund_groupedByReplicates.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}


##########################################################################
#
#   Copy all filtering with the collapsed asv count table
#
##########################################################################
asv_collapsedOn_taxonomy <- read.delim(as.character(args[9]), header=TRUE, stringsAsFactors=FALSE)

#Control filtering
controlFlag <- as.logical(args[6])
asv_count <- asv_collapsedOn_taxonomy

if (controlFlag == TRUE) {
  controls_df <- sample_metadata %>% select(Sample, controls)
  controls_df <- controls_df %>% filter(!is.na(controls))
  control_samples <- as.character(c(controls_df[1:nrow(controls_df), 1]))
  
  #TODO: CONTROL FILTERING (TBD)

  asv_count <- asv_count %>% select(-all_of(control_samples))
  asv_count$readsum <- rowSums(asv_count[,-1])
  asv_count <- asv_count %>% filter(readsum > 0)
  asv_count <- asv_count %>% select(-readsum)
  
  write.table(asv_count, "ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_controlsRemoved.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}

#Filter underperforming samples
filterLowQualSamplesFlag <- as.logical(args[7])
sampleFilterPerc <- as.numeric(args[8])

asv_count_samp <- asv_count
row.names(asv_count_samp) <- asv_count_samp$x
asv_count_samp <- asv_count_samp %>% select(-x)
asv_count_samp <- as.data.frame(t(asv_count_samp))
asv_count_samp <- tibble::rownames_to_column(asv_count_samp, "x")
asv_count_samp$readSum <- rowSums(asv_count_samp[,-1])
sample_sums <- as.numeric(asv_count_samp$readSum)
median_sample_effort <- median(sample_sums)

if (filterLowQualSamplesFlag == TRUE) {
  asv_count_samp <- asv_count_samp %>% filter(readSum >= median_sample_effort * (sampleFilterPerc/100))
  asv_count_samp <- asv_count_samp %>% select(-readSum)
  row.names(asv_count_samp) <- asv_count_samp$x
  asv_count_samp <- asv_count_samp %>% select(-x)
  asv_count_samp <- as.data.frame(t(asv_count_samp))
  asv_count_samp <- tibble::rownames_to_column(asv_count_samp, "x")
  asv_count <- asv_count_samp
  
  write.table(asv_count, "ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_lowEffortSamplesRemoved.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}

#Convert reads to rel. abund. (%)
asv_countrel <- asv_count
row.names(asv_countrel) <- asv_countrel$x
asv_countrel <- asv_countrel %>% select(-x)
asv_countrel <- as.data.frame(t(asv_countrel))
asv_countrel <- tibble::rownames_to_column(asv_countrel, "x")
asv_countrel$readSum <- rowSums(asv_countrel[,-1])
asv_countrel <- asv_countrel %>% mutate_at(vars(-x, -readSum), list(~(100*./readSum)))
asv_countrel <- asv_countrel %>% select(-readSum)
asv_countrel_rotated1 <- asv_countrel
asv_countrel_rotated2 <- asv_countrel
row.names(asv_countrel) <- asv_countrel$x
asv_countrel <- asv_countrel %>% select(-x)
asv_countrel <- as.data.frame(t(asv_countrel))
asv_countrel <- tibble::rownames_to_column(asv_countrel, "x")

write.table(asv_countrel, "ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)


#Grouped counts by mean (not sum) of rel. abund
replicateFlag <- as.logical(args[11])
sitelabelFlag <- as.logical(args[10])


if (sitelabelFlag == TRUE) { #Will expect "sites" column
  asv_countrel_rotated <- asv_countrel_rotated1
  asv_countrel_rotated <- asv_countrel_rotated %>% rename(Sample = x)
  asv_countrel_join <- left_join(select(sample_metadata, Sample, sites), asv_countrel_rotated, by = "Sample")
  temp_colnames <- colnames(asv_countrel_join)
  asv_countrel_join <- asv_countrel_join %>% filter(!is.na(eval(parse(text = temp_colnames[3]))))
  asv_countrel_join <- asv_countrel_join %>% select(-Sample)
  rep_group <- asv_countrel_join %>% group_by(sites) %>% 
    summarise_all(mean)
  rep_group <- rep_group %>% rename(x = sites)
  rep_group <- as.data.frame(rep_group)
  row.names(rep_group) <- rep_group$x
  rep_group <- rep_group %>% select(-x)
  rep_group <- as.data.frame(t(rep_group))
  rep_group <- tibble::rownames_to_column(rep_group, "x")
  write.table(rep_group, "ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund_groupedBySites.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}

if (replicateFlag == TRUE) { #Will expect "replicates" column
  asv_countrel_rotated <- asv_countrel_rotated2
  asv_countrel_rotated <- asv_countrel_rotated %>% rename(Sample = x)
  asv_countrel_join <- left_join(select(sample_metadata, Sample, replicates), asv_countrel_rotated, by = "Sample")
  temp_colnames <- colnames(asv_countrel_join)
  asv_countrel_join <- asv_countrel_join %>% filter(!is.na(eval(parse(text = temp_colnames[3]))))
  asv_countrel_join <- asv_countrel_join %>% select(-Sample)
  rep_group <- asv_countrel_join %>% group_by(replicates) %>% 
    summarise_all(mean)
  rep_group <- rep_group %>% rename(x = replicates)
  rep_group <- as.data.frame(rep_group)
  row.names(rep_group) <- rep_group$x
  rep_group <- rep_group %>% select(-x)
  rep_group <- as.data.frame(t(rep_group))
  rep_group <- tibble::rownames_to_column(rep_group, "x")
  write.table(rep_group, "ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund_groupedByReplicates.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}
