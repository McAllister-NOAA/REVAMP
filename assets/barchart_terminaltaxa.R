#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

########################################
#TEMP WHILE WORKING ON SCRIPT
# args[1]<-"/Users/mcallister/Desktop/test_figs" #FIGURE OUT directory
# args[2]<-"/Users/mcallister/Desktop/chris_test/18S/CP_all_out/ASV2Taxonomy/CP_all_out_barchart_forR.txt" #Location of barchart file from asv_taxonomy_processing_figureOuts.pl
# args[3]<-"/Users/mcallister/Desktop/chris_test/18S/CP_all_out/sample_order.txt" #sample order file
# args[4]<-"/Users/mcallister/Desktop/chris_test/18S/CP_all_out/ASV2Taxonomy/CP_all_out_barchart_forR_filtLowAbund_zzOther.txt" #filt low abund to zzOther file
########################################
library("ggplot2")
library("dplyr")
library("phyloseq")
library("ggpubr")

theme_set(theme_bw())

bar_chart <- read.delim(as.character(args[2]), header=TRUE, stringsAsFactors=FALSE)
sample_order <- read.delim(as.character(args[3]), header=FALSE, stringsAsFactors=FALSE)
sample_order <- as.vector(sample_order)

base_term_count <- as.numeric(system(paste0("cat ",as.character(args[2])," | cut -f3 | grep -v ","\"","TerminalTaxa","\""," | sort | uniq | wc -l", sep = ""),
                                       intern = TRUE))
w_choice <- 6.4 + ((base_term_count*2.3)/20)



setwd(paste0(as.character(args[1]), "/02_Barcharts/relative_abundance"))
bar_plot <- ggplot(bar_chart, aes(x = Sample, y = Value, fill = TerminalTaxa)) +
  geom_bar(position="fill", stat="identity") + 
  xlab("Samples") + ylab("Relative abundance") + labs(fill = "Unique Terminal Taxa") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid.major.x = element_blank())
bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = sample_order$V1)
bar_plot_legend <- get_legend(bar_plot)
bar_plot <- bar_plot + theme(legend.position='none')

pdf(file='barplot_relabund_allsamples_uniqueTerminalTaxa_withUnknowns_legend.pdf', width = w_choice, height = 6)
print(as_ggplot(bar_plot_legend))
dev.off()
pdf(file='barplot_relabund_allsamples_uniqueTerminalTaxa_withUnknowns.pdf', width = 11, height = 8.5)
print(bar_plot)
dev.off()

setwd(paste0(as.character(args[1]), "/02_Barcharts/read_count"))
bar_plot <- ggplot(bar_chart, aes(x = Sample, y = Value, fill = TerminalTaxa)) +
  geom_bar(position="stack", stat="identity") + 
  xlab("Samples") + ylab("Read Count") + labs(fill = "Unique Terminal Taxa") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid.major.x = element_blank())
bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = sample_order$V1)
bar_plot_legend <- get_legend(bar_plot)
bar_plot <- bar_plot + theme(legend.position='none')

pdf(file='barplot_readcount_allsamples_uniqueTerminalTaxa_withUnknowns_legend.pdf', width = w_choice, height = 6)
print(as_ggplot(bar_plot_legend))
dev.off()
pdf(file='barplot_readcount_allsamples_uniqueTerminalTaxa_withUnknowns.pdf', width = 11, height = 8.5)
print(bar_plot)
dev.off()


#Remove unknown taxa hits
bar_chart <- bar_chart %>% filter(TerminalTaxa != "Unknown")
bar_chart <- bar_chart %>% filter(TerminalTaxa != "Environmental Unknown")

setwd(paste0(as.character(args[1]), "/02_Barcharts/relative_abundance"))
bar_plot <- ggplot(bar_chart, aes(x = Sample, y = Value, fill = TerminalTaxa)) +
  geom_bar(position="fill", stat="identity") + 
  xlab("Samples") + ylab("Relative abundance") + labs(fill = "Unique Terminal Taxa") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid.major.x = element_blank())
bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = sample_order$V1)
bar_plot_legend <- get_legend(bar_plot)
bar_plot <- bar_plot + theme(legend.position='none')

pdf(file='barplot_relabund_allsamples_uniqueTerminalTaxa_noUnknowns_legend.pdf', width = w_choice, height = 6)
print(as_ggplot(bar_plot_legend))
dev.off()
pdf(file='barplot_relabund_allsamples_uniqueTerminalTaxa_noUnknowns.pdf', width = 11, height = 8.5)
print(bar_plot)
dev.off()

setwd(paste0(as.character(args[1]), "/02_Barcharts/read_count"))
bar_plot <- ggplot(bar_chart, aes(x = Sample, y = Value, fill = TerminalTaxa)) +
  geom_bar(position="stack", stat="identity") + 
  xlab("Samples") + ylab("Read Count") + labs(fill = "Unique Terminal Taxa") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid.major.x = element_blank())
bar_plot$data$Sample <- factor(bar_plot$data$Sample, levels = sample_order$V1)
bar_plot_legend <- get_legend(bar_plot)
bar_plot <- bar_plot + theme(legend.position='none')

pdf(file='barplot_readcount_allsamples_uniqueTerminalTaxa_noUnknowns_legend.pdf', width = w_choice, height = 6)
print(as_ggplot(bar_plot_legend))
dev.off()
pdf(file='barplot_readcount_allsamples_uniqueTerminalTaxa_noUnknowns.pdf', width = 11, height = 8.5)
print(bar_plot)
dev.off()

#Bar chart with filter low abundance taxa to zzOther
setwd(paste0(as.character(args[1]), "/02_Barcharts/relative_abundance"))
filter_bar_chart <- read.delim(as.character(args[4]), header=TRUE, stringsAsFactors=FALSE)

gg_color <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

count <- length(unique(filter_bar_chart$TerminalTaxa))
temp_col_pallete <- c(gg_color(count - 1), "grey")

filterbar_plot <- ggplot(filter_bar_chart, aes(x = Sample, y = Percent, fill = TerminalTaxa)) +
  geom_bar(position="stack", stat="identity") + 
  xlab("Samples") + ylab("Relative Abundance (%)") + labs(fill = "Unique Terminal Taxa") +
  scale_fill_manual(values = temp_col_pallete) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid.major.x = element_blank())
filterbar_plot$data$Sample <- factor(filterbar_plot$data$Sample, levels = sample_order$V1)
plot_legend <- get_legend(filterbar_plot)
filterbar_plot <- filterbar_plot + theme(legend.position='none')

pdf(file='barplot_relabund_allsamples_filtLowAbundTaxa_to_zzOther_uniqueTerminalTaxa_noUnknowns_legend.pdf', width = 11, height = 6)
print(as_ggplot(plot_legend))
dev.off()
pdf(file='barplot_relabund_allsamples_filtLowAbundTaxa_to_zzOther_uniqueTerminalTaxa_noUnknowns.pdf', width = 11, height = 8.5)
print(filterbar_plot)
dev.off()
