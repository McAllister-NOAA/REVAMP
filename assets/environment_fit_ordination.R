#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

########################################
#TEMP WHILE WORKING ON SCRIPT
# args[1]<-"/Users/mcallister/Desktop/test_figs" #FIGURE OUT directory
# args[2]<-"/Users/mcallister/Desktop/chris_test/18S/CP_all_out/processed_tables/ASVs_counts_NOUNKNOWNS_percentabund.tsv" #ASV rel abund
# args[3]<-"/Users/mcallister/Desktop/chris_test/18S/CP_all_out/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.tsv" #Taxonomy rel abund
# args[4]<-"/Users/mcallister/Desktop/chris_test/18S/CP_all_out/sample_metadata_forR.txt" #sample metadata
# args[5]<-"FALSE" #replicates flag
# args[6]<-"TRUE" #sites flag
# args[7]<-"TRUE" #chem data
# args[8]<-"/Users/mcallister/Desktop/chris_test/18S/CP_all_out/chem_headers.txt" #location chem headers
########################################
library("ggplot2")
library("dplyr")
library("vegan")
library("ggrepel")
library("ggalt")

setwd(as.character(args[1]))

theme_set(theme_bw())

replicateFlag <- as.logical(args[5])
sitelabelFlag <- as.logical(args[6])
chemDataFlag <- as.logical(args[7])


if (chemDataFlag == TRUE) {
  tryCatch({
  chem_headers <- read.delim(as.character(args[8]), header=FALSE, stringsAsFactors=FALSE)
  chem_headers <- as.vector(chem_headers)
  sample_metadata <- read.delim(as.character(args[4]), header=TRUE, stringsAsFactors=FALSE)
  reference_sample_metadata <- sample_metadata
  reference_sample_metadata2 <- sample_metadata
  sample_metadata <- sample_metadata %>% select(Sample, eval(parse(text = chem_headers)))
  
  ASV_percabund <- read.delim(as.character(args[2]), header=TRUE, stringsAsFactors=FALSE)
  row.names(ASV_percabund) <- ASV_percabund$x 
  ASV_percabund <- ASV_percabund %>% select(-x)
  ASV_percabund <- as.data.frame(t(ASV_percabund))
  
  ASV_meta <- ASV_percabund
  ASV_meta <- tibble::rownames_to_column(ASV_meta, "Sample")
  ASV_sample_meta <- right_join(sample_metadata, ASV_meta, by = "Sample")
  for (i in colnames(ASV_sample_meta)) {
    ASV_sample_meta <- ASV_sample_meta %>% filter(!is.na(eval(parse(text = i))))
    if (i != "Sample") {
      ASV_sample_meta <- ASV_sample_meta %>% filter(!is.na(as.numeric(eval(parse(text = i)))))
    }
  }
  row.names(ASV_sample_meta) <- ASV_sample_meta$Sample
  ASV_sample_meta <- ASV_sample_meta %>% select(-Sample)
  
  TAXA_percabund <- read.delim(as.character(args[3]), header=TRUE, stringsAsFactors=FALSE)
  row.names(TAXA_percabund) <- TAXA_percabund$x 
  TAXA_percabund <- TAXA_percabund %>% select(-x)
  TAXA_percabund <- as.data.frame(t(TAXA_percabund))
  
  TAXA_meta <- TAXA_percabund
  TAXA_meta <- tibble::rownames_to_column(TAXA_meta, "Sample")
  TAXA_sample_meta <- right_join(sample_metadata, TAXA_meta, by = "Sample")
  for (i in colnames(TAXA_sample_meta)) {
    TAXA_sample_meta <- TAXA_sample_meta %>% filter(!is.na(eval(parse(text = i))))
    if (i != "Sample") {
      TAXA_sample_meta <- TAXA_sample_meta %>% filter(!is.na(as.numeric(eval(parse(text = i)))))
    }
  }
  row.names(TAXA_sample_meta) <- TAXA_sample_meta$Sample
  TAXA_sample_meta <- TAXA_sample_meta %>% select(-Sample)
  
  # CHEMONLY_sample_meta <- sample_metadata
  # for (i in colnames(CHEMONLY_sample_meta)) {
  #   CHEMONLY_sample_meta <- CHEMONLY_sample_meta %>% filter(!is.na(eval(parse(text = i))))
  #   if (i != "Sample") {
  #     CHEMONLY_sample_meta <- CHEMONLY_sample_meta %>% filter(!is.na(as.numeric(eval(parse(text = i)))))
  #   }
  # }
  # row.names(CHEMONLY_sample_meta) <- CHEMONLY_sample_meta$Sample
  # CHEMONLY_sample_meta <- CHEMONLY_sample_meta %>% select(-Sample)
  
  }, error=function(e){})
} else {
  tryCatch({
  ASV_percabund <- read.delim(as.character(args[2]), header=TRUE, stringsAsFactors=FALSE)
  row.names(ASV_percabund) <- ASV_percabund$x 
  ASV_percabund <- ASV_percabund %>% select(-x)
  ASV_percabund <- as.data.frame(t(ASV_percabund))
  ASV_sample_meta <- ASV_percabund
  ASV_sample_meta <- tibble::rownames_to_column(ASV_sample_meta, "Sample")
  for (i in colnames(ASV_sample_meta)) {
    ASV_sample_meta <- ASV_sample_meta %>% filter(!is.na(eval(parse(text = i))))
    if (i != "Sample") {
      ASV_sample_meta <- ASV_sample_meta %>% filter(!is.na(as.numeric(eval(parse(text = i)))))
    }
  }
  row.names(ASV_sample_meta) <- ASV_sample_meta$Sample
  ASV_sample_meta <- ASV_sample_meta %>% select(-Sample)
  
  TAXA_percabund <- read.delim(as.character(args[3]), header=TRUE, stringsAsFactors=FALSE)
  row.names(TAXA_percabund) <- TAXA_percabund$x 
  TAXA_percabund <- TAXA_percabund %>% select(-x)
  TAXA_percabund <- as.data.frame(t(TAXA_percabund))
  TAXA_sample_meta <- TAXA_percabund
  TAXA_sample_meta <- tibble::rownames_to_column(TAXA_sample_meta, "Sample")
  for (i in colnames(TAXA_sample_meta)) {
    TAXA_sample_meta <- TAXA_sample_meta %>% filter(!is.na(eval(parse(text = i))))
    if (i != "Sample") {
      TAXA_sample_meta <- TAXA_sample_meta %>% filter(!is.na(as.numeric(eval(parse(text = i)))))
    }
  }
  row.names(TAXA_sample_meta) <- TAXA_sample_meta$Sample
  TAXA_sample_meta <- TAXA_sample_meta %>% select(-Sample)
  
  sample_metadata <- read.delim(as.character(args[4]), header=TRUE, stringsAsFactors=FALSE)
  reference_sample_metadata <- sample_metadata
  reference_sample_metadata2 <- sample_metadata
  }, error=function(e){})
}

#ASV-based ordination
setwd(paste0(as.character(args[1]),"/08_EnvironmentFit_Ordination/ASV_based"))
tryCatch({
NMDS_ASV <- metaMDS(ASV_percabund, k=2, trymax = 100, trace = F, autotransform = FALSE, distance = "bray")
NMDS_tab <- as.data.frame(NMDS_ASV$species)
NMDS_tab <- tibble::rownames_to_column(NMDS_tab, "ASV")
write.table(NMDS_tab, "NMDS_vegan_ASVbased_ASV_coordinates.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
NMDS_tab <- as.data.frame(NMDS_ASV$points)
NMDS_tab <- tibble::rownames_to_column(NMDS_tab, "Sample")
write.table(NMDS_tab, "NMDS_vegan_ASVbased_sample_coordinates.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
reference_sample_metadata <- left_join(reference_sample_metadata, NMDS_tab, by = "Sample")
reference_sample_metadata <- reference_sample_metadata %>% filter(!is.na(MDS1))
}, error=function(e){})

tryCatch({
pdf(file='NMDS_vegan_ASVbased_relabund_stressplot.pdf', width = 11, height = 8.5)
stressplot(NMDS_ASV)
dev.off()
}, error=function(e){})

tryCatch({
fitenv_ASV <- envfit(NMDS_ASV, ASV_sample_meta, permutations = 999)
arrows <- as.data.frame(fitenv_ASV$vectors$arrows*sqrt(fitenv_ASV$vectors$r))
arrows <- tibble::rownames_to_column(arrows, "Variable")
r2 <- as.data.frame(fitenv_ASV$vectors$r)
r2 <- r2 %>% rename(r_square = "fitenv_ASV$vectors$r")
r2 <- tibble::rownames_to_column(r2, "Variable")
pvals <- as.data.frame(fitenv_ASV$vectors$pvals)
pvals <- pvals %>% rename(pval = "fitenv_ASV$vectors$pvals")
pvals <- tibble::rownames_to_column(pvals, "Variable")
filtenv_tab <- left_join(arrows, r2, by = "Variable")
filtenv_tab <- left_join(filtenv_tab, pvals, by = "Variable")
write.table(filtenv_tab, "NMDS_vegan_ASVbased_environmentFitVectors_coordinates_r2_pvalues.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}, error=function(e){})

tryCatch({
ord_plot <- ggplot(NMDS_tab) + geom_point(aes(x = MDS1, y = MDS2), color="red") +
  geom_text_repel(aes(label = Sample, x = MDS1, y = MDS2), box.padding = 0.4, min.segment.length = 0.2) +
  xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank())
pdf(file='NMDS_vegan_ASVbased_samples_relabund_filtsamples_alltaxa.pdf', width = 11, height = 8.5)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
filtenv_plot <- filtenv_tab %>% filter(pval <= 0.05)
ord_plot <- ggplot(NMDS_tab) + geom_point(aes(x = MDS1, y = MDS2), color="red") +
  geom_segment(data = filtenv_plot, aes(x=0,xend= NMDS1, y=0, yend = NMDS2), 
               arrow = arrow(length=unit(0.2, "cm")), 
               color="royalblue2", 
               inherit.aes = FALSE) +
  geom_text(data = filtenv_plot, aes(x=NMDS1, y=NMDS2, label=Variable), color="royalblue4") +
  xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank())
pdf(file='NMDS_vegan_ASVbased_samples_relabund_filtsamples_alltaxa_environmentFitVectors_p0.05.pdf', width = 11, height = 8.5)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
filtenv_plot <- filtenv_tab %>% filter(pval <= 0.001)
ord_plot <- ggplot(NMDS_tab) + geom_point(aes(x = MDS1, y = MDS2), color="red") +
  geom_segment(data = filtenv_plot, aes(x=0,xend= NMDS1, y=0, yend = NMDS2), 
               arrow = arrow(length=unit(0.2, "cm")), 
               color="royalblue2", 
               inherit.aes = FALSE) +
  geom_text(data = filtenv_plot, aes(x=NMDS1, y=NMDS2, label=Variable), color="royalblue4") +
  xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank())
pdf(file='NMDS_vegan_ASVbased_samples_relabund_filtsamples_alltaxa_environmentFitVectors_p0.001.pdf', width = 11, height = 8.5)
print(ord_plot)
dev.off()
}, error=function(e){})

if (chemDataFlag == TRUE) {
  tryCatch({
filtenv_plot <- filtenv_tab %>% filter(Variable %in% chem_headers$V1)
ord_plot <- ggplot(NMDS_tab) + geom_point(aes(x = MDS1, y = MDS2), color="red") +
  geom_segment(data = filtenv_plot, aes(x=0,xend= NMDS1, y=0, yend = NMDS2), 
               arrow = arrow(length=unit(0.2, "cm")), 
               color="royalblue2", 
               inherit.aes = FALSE) +
  geom_text(data = filtenv_plot, aes(x=NMDS1, y=NMDS2, label=Variable), color="royalblue4") +
  xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank())
pdf(file='NMDS_vegan_ASVbased_samples_relabund_filtsamples_alltaxa_environmentFitVectors_chemOnly.pdf', width = 11, height = 8.5)
print(ord_plot)
dev.off()
}, error=function(e){})
}

if (replicateFlag == TRUE) {
  tryCatch({
  ord_plot <- ggplot(reference_sample_metadata) + geom_point(aes(x = MDS1, y = MDS2, color = replicates)) +
    xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank()) + 
    geom_encircle(aes(x = MDS1, y = MDS2, fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.3, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='NMDS_vegan_ASVbased_samples_relabund_filtsamples_alltaxa_encircleReplicates.pdf', width = 11, height = 8.5)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  filtenv_plot <- filtenv_tab %>% filter(pval <= 0.05)
  ord_plot <- ggplot(reference_sample_metadata) + geom_point(aes(x = MDS1, y = MDS2, color = replicates)) +
    geom_encircle(aes(x = MDS1, y = MDS2, fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.3, expand = 0, s_shape = 1, spread = 0.001) +
    geom_segment(data = filtenv_plot, aes(x=0,xend= NMDS1, y=0, yend = NMDS2), 
                 arrow = arrow(length=unit(0.2, "cm")), 
                 color="royalblue2", 
                 inherit.aes = FALSE) +
    geom_text(data = filtenv_plot, aes(x=NMDS1, y=NMDS2, label=Variable), color="royalblue4") +
    xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank())
  pdf(file='NMDS_vegan_ASVbased_samples_relabund_filtsamples_alltaxa_environmentFitVectors_p0.05_encircleReplicates.pdf', width = 11, height = 8.5)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  filtenv_plot <- filtenv_tab %>% filter(pval <= 0.001)
  ord_plot <- ggplot(reference_sample_metadata) + geom_point(aes(x = MDS1, y = MDS2, color = replicates)) +
    geom_encircle(aes(x = MDS1, y = MDS2, fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.3, expand = 0, s_shape = 1, spread = 0.001) +
    geom_segment(data = filtenv_plot, aes(x=0,xend= NMDS1, y=0, yend = NMDS2), 
                 arrow = arrow(length=unit(0.2, "cm")), 
                 color="royalblue2", 
                 inherit.aes = FALSE) +
    geom_text(data = filtenv_plot, aes(x=NMDS1, y=NMDS2, label=Variable), color="royalblue4") +
    xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank())
  pdf(file='NMDS_vegan_ASVbased_samples_relabund_filtsamples_alltaxa_environmentFitVectors_p0.001_encircleReplicates.pdf', width = 11, height = 8.5)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  if (chemDataFlag == TRUE) {
    tryCatch({
    filtenv_plot <- filtenv_tab %>% filter(Variable %in% chem_headers$V1)
    ord_plot <- ggplot(reference_sample_metadata) + geom_point(aes(x = MDS1, y = MDS2, color = replicates)) +
      geom_encircle(aes(x = MDS1, y = MDS2, fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.3, expand = 0, s_shape = 1, spread = 0.001) +
      geom_segment(data = filtenv_plot, aes(x=0,xend= NMDS1, y=0, yend = NMDS2), 
                   arrow = arrow(length=unit(0.2, "cm")), 
                   color="royalblue2", 
                   inherit.aes = FALSE) +
      geom_text(data = filtenv_plot, aes(x=NMDS1, y=NMDS2, label=Variable), color="royalblue4") +
      xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank())
    pdf(file='NMDS_vegan_ASVbased_samples_relabund_filtsamples_alltaxa_environmentFitVectors_chemOnly_encircleReplicates.pdf', width = 11, height = 8.5)
    print(ord_plot)
    dev.off()
    }, error=function(e){})
  }
}

if (sitelabelFlag == TRUE) {
  tryCatch({
  ord_plot <- ggplot(reference_sample_metadata) + geom_point(aes(x = MDS1, y = MDS2, color = sites)) +
    xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank()) + 
    geom_encircle(aes(x = MDS1, y = MDS2, fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.3, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='NMDS_vegan_ASVbased_samples_relabund_filtsamples_alltaxa_encircleSites.pdf', width = 11, height = 8.5)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  filtenv_plot <- filtenv_tab %>% filter(pval <= 0.05)
  ord_plot <- ggplot(reference_sample_metadata) + geom_point(aes(x = MDS1, y = MDS2, color = sites)) +
    geom_encircle(aes(x = MDS1, y = MDS2, fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.3, expand = 0, s_shape = 1, spread = 0.001) +
    geom_segment(data = filtenv_plot, aes(x=0,xend= NMDS1, y=0, yend = NMDS2), 
                 arrow = arrow(length=unit(0.2, "cm")), 
                 color="royalblue2", 
                 inherit.aes = FALSE) +
    geom_text(data = filtenv_plot, aes(x=NMDS1, y=NMDS2, label=Variable), color="royalblue4") +
    xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank())
  pdf(file='NMDS_vegan_ASVbased_samples_relabund_filtsamples_alltaxa_environmentFitVectors_p0.05_encircleSites.pdf', width = 11, height = 8.5)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  filtenv_plot <- filtenv_tab %>% filter(pval <= 0.001)
  ord_plot <- ggplot(reference_sample_metadata) + geom_point(aes(x = MDS1, y = MDS2, color = sites)) +
    geom_encircle(aes(x = MDS1, y = MDS2, fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.3, expand = 0, s_shape = 1, spread = 0.001) +
    geom_segment(data = filtenv_plot, aes(x=0,xend= NMDS1, y=0, yend = NMDS2), 
                 arrow = arrow(length=unit(0.2, "cm")), 
                 color="royalblue2", 
                 inherit.aes = FALSE) +
    geom_text(data = filtenv_plot, aes(x=NMDS1, y=NMDS2, label=Variable), color="royalblue4") +
    xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank())
  pdf(file='NMDS_vegan_ASVbased_samples_relabund_filtsamples_alltaxa_environmentFitVectors_p0.001_encircleSites.pdf', width = 11, height = 8.5)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  if (chemDataFlag == TRUE) {
    tryCatch({
    filtenv_plot <- filtenv_tab %>% filter(Variable %in% chem_headers$V1)
    ord_plot <- ggplot(reference_sample_metadata) + geom_point(aes(x = MDS1, y = MDS2, color = sites)) +
      geom_encircle(aes(x = MDS1, y = MDS2, fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.3, expand = 0, s_shape = 1, spread = 0.001) +
      geom_segment(data = filtenv_plot, aes(x=0,xend= NMDS1, y=0, yend = NMDS2), 
                   arrow = arrow(length=unit(0.2, "cm")), 
                   color="royalblue2", 
                   inherit.aes = FALSE) +
      geom_text(data = filtenv_plot, aes(x=NMDS1, y=NMDS2, label=Variable), color="royalblue4") +
      xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank())
    pdf(file='NMDS_vegan_ASVbased_samples_relabund_filtsamples_alltaxa_environmentFitVectors_chemOnly_encircleSites.pdf', width = 11, height = 8.5)
    print(ord_plot)
    dev.off()
    }, error=function(e){})
  }
}

#Taxonomy-based ordination
setwd(paste0(as.character(args[1]),"/08_EnvironmentFit_Ordination/Taxonomy_merge_based"))

tryCatch({
NMDS_TAXA <- metaMDS(TAXA_percabund, k=2, trymax = 100, trace = F, autotransform = FALSE, distance = "bray")
NMDS_tab <- as.data.frame(NMDS_TAXA$species)
NMDS_tab <- tibble::rownames_to_column(NMDS_tab, "ASV")
write.table(NMDS_tab, "NMDS_vegan_TAXAbased_ASV_coordinates.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
NMDS_tab <- as.data.frame(NMDS_TAXA$points)
NMDS_tab <- tibble::rownames_to_column(NMDS_tab, "Sample")
write.table(NMDS_tab, "NMDS_vegan_TAXAbased_sample_coordinates.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}, error=function(e){})

tryCatch({
reference_sample_metadata <- reference_sample_metadata2
reference_sample_metadata <- left_join(reference_sample_metadata, NMDS_tab, by = "Sample")
reference_sample_metadata <- reference_sample_metadata %>% filter(!is.na(MDS1))
}, error=function(e){})

tryCatch({
pdf(file='NMDS_vegan_TAXAbased_relabund_stressplot.pdf', width = 11, height = 8.5)
stressplot(NMDS_TAXA)
dev.off()
}, error=function(e){})

tryCatch({
fitenv_TAXA <- envfit(NMDS_TAXA, TAXA_sample_meta, permutations = 999)
arrows <- as.data.frame(fitenv_TAXA$vectors$arrows*sqrt(fitenv_TAXA$vectors$r))
arrows <- tibble::rownames_to_column(arrows, "Variable")
r2 <- as.data.frame(fitenv_TAXA$vectors$r)
r2 <- r2 %>% rename(r_square = "fitenv_TAXA$vectors$r")
r2 <- tibble::rownames_to_column(r2, "Variable")
pvals <- as.data.frame(fitenv_TAXA$vectors$pvals)
pvals <- pvals %>% rename(pval = "fitenv_TAXA$vectors$pvals")
pvals <- tibble::rownames_to_column(pvals, "Variable")
filtenv_tab <- left_join(arrows, r2, by = "Variable")
filtenv_tab <- left_join(filtenv_tab, pvals, by = "Variable")
write.table(filtenv_tab, "NMDS_vegan_TAXAbased_environmentFitVectors_coordinates_r2_pvalues.tsv", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}, error=function(e){})

tryCatch({
ord_plot <- ggplot(NMDS_tab) + geom_point(aes(x = MDS1, y = MDS2), color="red") +
  geom_text_repel(aes(label = Sample, x = MDS1, y = MDS2), box.padding = 0.4, min.segment.length = 0.2) +
  xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank())
pdf(file='NMDS_vegan_TAXAbased_samples_relabund_filtsamples_alltaxa.pdf', width = 11, height = 8.5)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
filtenv_plot <- filtenv_tab %>% filter(pval <= 0.05)
ord_plot <- ggplot(NMDS_tab) + geom_point(aes(x = MDS1, y = MDS2), color="red") +
  geom_segment(data = filtenv_plot, aes(x=0,xend= NMDS1, y=0, yend = NMDS2), 
               arrow = arrow(length=unit(0.2, "cm")), 
               color="royalblue2", 
               inherit.aes = FALSE) +
  geom_text(data = filtenv_plot, aes(x=NMDS1, y=NMDS2, label=Variable), color="royalblue4") +
  xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank())
pdf(file='NMDS_vegan_TAXAbased_samples_relabund_filtsamples_alltaxa_environmentFitVectors_p0.05.pdf', width = 11, height = 8.5)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
filtenv_plot <- filtenv_tab %>% filter(pval <= 0.001)
ord_plot <- ggplot(NMDS_tab) + geom_point(aes(x = MDS1, y = MDS2), color="red") +
  geom_segment(data = filtenv_plot, aes(x=0,xend= NMDS1, y=0, yend = NMDS2), 
               arrow = arrow(length=unit(0.2, "cm")), 
               color="royalblue2", 
               inherit.aes = FALSE) +
  geom_text(data = filtenv_plot, aes(x=NMDS1, y=NMDS2, label=Variable), color="royalblue4") +
  xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank())
pdf(file='NMDS_vegan_TAXAbased_samples_relabund_filtsamples_alltaxa_environmentFitVectors_p0.001.pdf', width = 11, height = 8.5)
print(ord_plot)
dev.off()
}, error=function(e){})

if (chemDataFlag == TRUE) {
  tryCatch({
  filtenv_plot <- filtenv_tab %>% filter(Variable %in% chem_headers$V1)
  ord_plot <- ggplot(NMDS_tab) + geom_point(aes(x = MDS1, y = MDS2), color="red") +
    geom_segment(data = filtenv_plot, aes(x=0,xend= NMDS1, y=0, yend = NMDS2), 
                 arrow = arrow(length=unit(0.2, "cm")), 
                 color="royalblue2", 
                 inherit.aes = FALSE) +
    geom_text(data = filtenv_plot, aes(x=NMDS1, y=NMDS2, label=Variable), color="royalblue4") +
    xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank())
  pdf(file='NMDS_vegan_TAXAbased_samples_relabund_filtsamples_alltaxa_environmentFitVectors_chemOnly.pdf', width = 11, height = 8.5)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
}


if (replicateFlag == TRUE) {
  tryCatch({
  ord_plot <- ggplot(reference_sample_metadata) + geom_point(aes(x = MDS1, y = MDS2, color = replicates)) +
    xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank()) + 
    geom_encircle(aes(x = MDS1, y = MDS2, fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.3, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='NMDS_vegan_TAXAbased_samples_relabund_filtsamples_alltaxa_encircleReplicates.pdf', width = 11, height = 8.5)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  filtenv_plot <- filtenv_tab %>% filter(pval <= 0.05)
  ord_plot <- ggplot(reference_sample_metadata) + geom_point(aes(x = MDS1, y = MDS2, color = replicates)) +
    geom_encircle(aes(x = MDS1, y = MDS2, fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.3, expand = 0, s_shape = 1, spread = 0.001) +
    geom_segment(data = filtenv_plot, aes(x=0,xend= NMDS1, y=0, yend = NMDS2), 
                 arrow = arrow(length=unit(0.2, "cm")), 
                 color="royalblue2", 
                 inherit.aes = FALSE) +
    geom_text(data = filtenv_plot, aes(x=NMDS1, y=NMDS2, label=Variable), color="royalblue4") +
    xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank())
  pdf(file='NMDS_vegan_TAXAbased_samples_relabund_filtsamples_alltaxa_environmentFitVectors_p0.05_encircleReplicates.pdf', width = 11, height = 8.5)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  filtenv_plot <- filtenv_tab %>% filter(pval <= 0.001)
  ord_plot <- ggplot(reference_sample_metadata) + geom_point(aes(x = MDS1, y = MDS2, color = replicates)) +
    geom_encircle(aes(x = MDS1, y = MDS2, fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.3, expand = 0, s_shape = 1, spread = 0.001) +
    geom_segment(data = filtenv_plot, aes(x=0,xend= NMDS1, y=0, yend = NMDS2), 
                 arrow = arrow(length=unit(0.2, "cm")), 
                 color="royalblue2", 
                 inherit.aes = FALSE) +
    geom_text(data = filtenv_plot, aes(x=NMDS1, y=NMDS2, label=Variable), color="royalblue4") +
    xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank())
  pdf(file='NMDS_vegan_TAXAbased_samples_relabund_filtsamples_alltaxa_environmentFitVectors_p0.001_encircleReplicates.pdf', width = 11, height = 8.5)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  if (chemDataFlag == TRUE) {
    tryCatch({
    filtenv_plot <- filtenv_tab %>% filter(Variable %in% chem_headers$V1)
    ord_plot <- ggplot(reference_sample_metadata) + geom_point(aes(x = MDS1, y = MDS2, color = replicates)) +
      geom_encircle(aes(x = MDS1, y = MDS2, fill = replicates, group = replicates), inherit.aes = TRUE, alpha = 0.3, expand = 0, s_shape = 1, spread = 0.001) +
      geom_segment(data = filtenv_plot, aes(x=0,xend= NMDS1, y=0, yend = NMDS2), 
                   arrow = arrow(length=unit(0.2, "cm")), 
                   color="royalblue2", 
                   inherit.aes = FALSE) +
      geom_text(data = filtenv_plot, aes(x=NMDS1, y=NMDS2, label=Variable), color="royalblue4") +
      xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank())
    pdf(file='NMDS_vegan_TAXAbased_samples_relabund_filtsamples_alltaxa_environmentFitVectors_chemOnly_encircleReplicates.pdf', width = 11, height = 8.5)
    print(ord_plot)
    dev.off()
    }, error=function(e){})
  }
}

if (sitelabelFlag == TRUE) {
  tryCatch({
  ord_plot <- ggplot(reference_sample_metadata) + geom_point(aes(x = MDS1, y = MDS2, color = sites)) +
    xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank()) + 
    geom_encircle(aes(x = MDS1, y = MDS2, fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.3, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='NMDS_vegan_TAXAbased_samples_relabund_filtsamples_alltaxa_encircleSites.pdf', width = 11, height = 8.5)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  filtenv_plot <- filtenv_tab %>% filter(pval <= 0.05)
  ord_plot <- ggplot(reference_sample_metadata) + geom_point(aes(x = MDS1, y = MDS2, color = sites)) +
    geom_encircle(aes(x = MDS1, y = MDS2, fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.3, expand = 0, s_shape = 1, spread = 0.001) +
    geom_segment(data = filtenv_plot, aes(x=0,xend= NMDS1, y=0, yend = NMDS2), 
                 arrow = arrow(length=unit(0.2, "cm")), 
                 color="royalblue2", 
                 inherit.aes = FALSE) +
    geom_text(data = filtenv_plot, aes(x=NMDS1, y=NMDS2, label=Variable), color="royalblue4") +
    xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank())
  pdf(file='NMDS_vegan_TAXAbased_samples_relabund_filtsamples_alltaxa_environmentFitVectors_p0.05_encircleSites.pdf', width = 11, height = 8.5)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  tryCatch({
  filtenv_plot <- filtenv_tab %>% filter(pval <= 0.001)
  ord_plot <- ggplot(reference_sample_metadata) + geom_point(aes(x = MDS1, y = MDS2, color = sites)) +
    geom_encircle(aes(x = MDS1, y = MDS2, fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.3, expand = 0, s_shape = 1, spread = 0.001) +
    geom_segment(data = filtenv_plot, aes(x=0,xend= NMDS1, y=0, yend = NMDS2), 
                 arrow = arrow(length=unit(0.2, "cm")), 
                 color="royalblue2", 
                 inherit.aes = FALSE) +
    geom_text(data = filtenv_plot, aes(x=NMDS1, y=NMDS2, label=Variable), color="royalblue4") +
    xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank())
  pdf(file='NMDS_vegan_TAXAbased_samples_relabund_filtsamples_alltaxa_environmentFitVectors_p0.001_encircleSites.pdf', width = 11, height = 8.5)
  print(ord_plot)
  dev.off()
  }, error=function(e){})
  
  if (chemDataFlag == TRUE) {
    tryCatch({
    filtenv_plot <- filtenv_tab %>% filter(Variable %in% chem_headers$V1)
    ord_plot <- ggplot(reference_sample_metadata) + geom_point(aes(x = MDS1, y = MDS2, color = sites)) +
      geom_encircle(aes(x = MDS1, y = MDS2, fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.3, expand = 0, s_shape = 1, spread = 0.001) +
      geom_segment(data = filtenv_plot, aes(x=0,xend= NMDS1, y=0, yend = NMDS2), 
                   arrow = arrow(length=unit(0.2, "cm")), 
                   color="royalblue2", 
                   inherit.aes = FALSE) +
      geom_text(data = filtenv_plot, aes(x=NMDS1, y=NMDS2, label=Variable), color="royalblue4") +
      xlab("NMDS1") + ylab("NMDS2") + theme(panel.grid = element_blank())
    pdf(file='NMDS_vegan_TAXAbased_samples_relabund_filtsamples_alltaxa_environmentFitVectors_chemOnly_encircleSites.pdf', width = 11, height = 8.5)
    print(ord_plot)
    dev.off()
    }, error=function(e){})
  }
}
