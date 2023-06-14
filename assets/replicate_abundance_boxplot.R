#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

########################################
# # TEMP WHILE WORKING ON SCRIPT
# args[1]<-"/PATH/REVAMP/outdir/test_figs" #FIGURE OUT directory
# args[2]<-"/PATH/REVAMP/outdir/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_ASVbased_withUnknowns_filtsamples.txt" #presenceabsence_compRelAbund_x file
# args[3]<-"/PATH/REVAMP/outdir/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_ASVbased_NoUnknowns_filtsamples.txt" #presenceabsence_compRelAbund_x file
# args[4]<-"/PATH/REVAMP/outdir/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_TAXAbased_NoUnknowns_filtsamples.txt" #presenceabsence_compRelAbund_x file

########################################
library("dplyr")
library("ggplot2")

setwd(as.character(args[1]))

theme_set(theme_bw())

PA_table <- read.delim(as.character(args[2]), header=TRUE, stringsAsFactors=FALSE) #ASV based, with unknowns

PA_table <- PA_table %>% filter(num_present > 0)
PA_table$num_present <- as.character(PA_table$num_present)
PA_table$sum_counts <- log(PA_table$sum_counts, 2)
PA_table$avg_counts_wherefound <- log(PA_table$avg_counts_wherefound, 2)
PA_table$avg_counts_allrep <- log(PA_table$avg_counts_allrep, 2)


PA_plot <- ggplot(PA_table, aes(num_present, sum_counts)) +
  geom_violin(scale = "count", fill = "plum") +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1) +
  labs(y = "log2(Total reads per ASV per site)",
       x = "Number of replicates ASV detected per site")

pdf(file='ViolinBoxPlot_ReadsVSReplicateDetection_ASVbased_withUnknowns_sumReads.pdf', width = 11, height = 8.5)
print(PA_plot)
dev.off()

PA_plot <- ggplot(PA_table, aes(num_present, avg_counts_wherefound)) +
  geom_violin(scale = "count", fill = "plum") +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1) +
  labs(y = "log2(Average reads over replicates detected per ASV per site)",
       x = "Number of replicates ASV detected per site")

pdf(file='ViolinBoxPlot_ReadsVSReplicateDetection_ASVbased_withUnknowns_avgReadsWhereDeteted.pdf', width = 11, height = 8.5)
print(PA_plot)
dev.off()

PA_plot <- ggplot(PA_table, aes(num_present, avg_counts_allrep)) +
  geom_violin(scale = "count", fill = "plum") +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1) +
  labs(y = "log2(Average reads over all replicates per ASV per site)",
       x = "Number of replicates ASV detected per site")

pdf(file='ViolinBoxPlot_ReadsVSReplicateDetection_ASVbased_withUnknowns_avgReadsAllReplicates.pdf', width = 11, height = 8.5)
print(PA_plot)
dev.off()


#############################################################################

PA_table <- read.delim(as.character(args[3]), header=TRUE, stringsAsFactors=FALSE) #ASV based, without unknowns

PA_table <- PA_table %>% filter(num_present > 0)
PA_table$num_present <- as.character(PA_table$num_present)
PA_table$sum_counts <- log(PA_table$sum_counts, 2)
PA_table$avg_counts_wherefound <- log(PA_table$avg_counts_wherefound, 2)
PA_table$avg_counts_allrep <- log(PA_table$avg_counts_allrep, 2)


PA_plot <- ggplot(PA_table, aes(num_present, sum_counts)) +
  geom_violin(scale = "count", fill = "plum") +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1) +
  labs(y = "log2(Total reads per ASV per site)",
       x = "Number of replicates ASV detected per site")

pdf(file='ViolinBoxPlot_ReadsVSReplicateDetection_ASVbased_withoutUnknowns_sumReads.pdf', width = 11, height = 8.5)
print(PA_plot)
dev.off()

PA_plot <- ggplot(PA_table, aes(num_present, avg_counts_wherefound)) +
  geom_violin(scale = "count", fill = "plum") +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1) +
  labs(y = "log2(Average reads over replicates detected per ASV per site)",
       x = "Number of replicates ASV detected per site")

pdf(file='ViolinBoxPlot_ReadsVSReplicateDetection_ASVbased_withoutUnknowns_avgReadsWhereDeteted.pdf', width = 11, height = 8.5)
print(PA_plot)
dev.off()

PA_plot <- ggplot(PA_table, aes(num_present, avg_counts_allrep)) +
  geom_violin(scale = "count", fill = "plum") +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1) +
  labs(y = "log2(Average reads over all replicates per ASV per site)",
       x = "Number of replicates ASV detected per site")

pdf(file='ViolinBoxPlot_ReadsVSReplicateDetection_ASVbased_withoutUnknowns_avgReadsAllReplicates.pdf', width = 11, height = 8.5)
print(PA_plot)
dev.off()

#############################################################################

PA_table <- read.delim(as.character(args[4]), header=TRUE, stringsAsFactors=FALSE) #TAXA based, without unknowns

PA_table <- PA_table %>% filter(num_present > 0)
PA_table$num_present <- as.character(PA_table$num_present)
PA_table$sum_counts <- log(PA_table$sum_counts, 2)
PA_table$avg_counts_wherefound <- log(PA_table$avg_counts_wherefound, 2)
PA_table$avg_counts_allrep <- log(PA_table$avg_counts_allrep, 2)


PA_plot <- ggplot(PA_table, aes(num_present, sum_counts)) +
  geom_violin(scale = "count", fill = "plum") +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1) +
  labs(y = "log2(Total reads per taxa per site)",
       x = "Number of replicates taxa detected per site")

pdf(file='ViolinBoxPlot_ReadsVSReplicateDetection_TAXAbased_withoutUnknowns_sumReads.pdf', width = 11, height = 8.5)
print(PA_plot)
dev.off()

PA_plot <- ggplot(PA_table, aes(num_present, avg_counts_wherefound)) +
  geom_violin(scale = "count", fill = "plum") +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1) +
  labs(y = "log2(Average reads over replicates detected per taxa per site)",
       x = "Number of replicates taxa detected per site")

pdf(file='ViolinBoxPlot_ReadsVSReplicateDetection_TAXAbased_withoutUnknowns_avgReadsWhereDeteted.pdf', width = 11, height = 8.5)
print(PA_plot)
dev.off()

PA_plot <- ggplot(PA_table, aes(num_present, avg_counts_allrep)) +
  geom_violin(scale = "count", fill = "plum") +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1) +
  labs(y = "log2(Average reads over all replicates per taxa per site)",
       x = "Number of replicates taxa detected per site")

pdf(file='ViolinBoxPlot_ReadsVSReplicateDetection_TAXAbased_withoutUnknowns_avgReadsAllReplicates.pdf', width = 11, height = 8.5)
print(PA_plot)
dev.off()
