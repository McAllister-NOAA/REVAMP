#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

########################################
#TEMP WHILE WORKING ON SCRIPT
# args[1]<-"" #FIGURE OUT directory
# args[2]<-"" #MergedMarkers_asvTaxonomyTable_NOUNKNOWNS.txt
# args[3]<-"" #MergedMarkers_ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.txt
# args[4]<-"" #MergedMarkers_sample_metadata_forR.txt
# args[5]<-"FALSE" #replicateFlag
# args[6]<-"TRUE" #sitelabelFlag
# args[7]<-"TRUE" #groups defined in sample metadata file
# args[8]<-"3" #number of groups defined in sample metadata file
# args[9]<-"TRUE" #whether taxa of interest file is given
# args[10]<-"Order" #category to filter on for taxa of interest
# args[11]<-"" #location of taxa of interest one per line
########################################
library("ggplot2")
library("dplyr")
library("phyloseq")
library("ggpubr")
library("ggalt")
library("ggrepel")
library("vegan")
library("spatstat")

setwd(as.character(args[1]))
theme_set(theme_bw())

replicateFlag <- as.logical(args[5])
sitelabelFlag <- as.logical(args[6])
taxaofinterestFlag <- as.logical(args[9])
groupingFlag <- as.logical(args[7])
numberofGroups <- as.numeric(args[8])

##################################
#
#  Import Taxonomy-based, relative abund, qual filtered (phylo11)
#
##################################
asv_count <- read.delim(as.character(args[3]), header=TRUE, stringsAsFactors=FALSE)
row.names(asv_count) <- asv_count$x
asv_count <- asv_count %>% select(-x)
asv_count_mat <- as.matrix(asv_count)

asv_taxonomy <- read.delim(as.character(args[2]), header=TRUE, stringsAsFactors=FALSE)
row.names(asv_taxonomy) <- asv_taxonomy$ASV
asv_taxonomy <- asv_taxonomy %>% select(-ASV)
asv_taxonomy_mat <- as.matrix(asv_taxonomy)

sample_metadata <- read.delim(as.character(args[4]), header=TRUE, stringsAsFactors=TRUE)
row.names(sample_metadata) <- sample_metadata$Sample
sample_metadata <- sample_metadata %>% select(-Sample)

ASV = otu_table(asv_count_mat, taxa_are_rows = TRUE)
TAX = tax_table(asv_taxonomy_mat)
samples = sample_data(sample_metadata)

phylo11 <- phyloseq(ASV, TAX, samples)

message("#######################################################################################")
message("#")
message("# Phyloseq objects have been created to work with:")
message("#")
message("# Taxonomy-based, relative abund, qual filtered (phylo11) ALL")
message("#")
message("#######################################################################################")

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

options(max.print=1000000)

#Create ordinations for phyloseq object (NMDS & PCoA)
tryCatch({
message("phylo11 NMDS")
setwd(paste0(as.character(args[1]), "/Figures/Ordination", sep = ""))
phylo11_NMDS.ord <- ordinate(phylo11, "NMDS")
phylo11_NMDS.ord
phylo11_stress <- phylo11_NMDS.ord$stress
ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, type="taxa")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
write.table(ord_tab, "NMDS_TAXAbased_taxa_relabund_allsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, type="samples")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "Samples")
ord_tab <- ord_tab %>% select(Samples, NMDS1, NMDS2)
write.table(ord_tab, "NMDS_TAXAbased_samples_relabund_allsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
phylo11_PCoA.ord <- ordinate(phylo11, "PCoA")
message("phylo11 PCoA")
phylo11_PCoA.ord
ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, type="taxa")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "ASV")
write.table(ord_tab, "PCoA_TAXAbased_taxa_relabund_allsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, type="samples")
ord_tab <- ord_plot$data
ord_tab <- tibble::rownames_to_column(ord_tab, "Samples")
ord_tab <- ord_tab %>% select(Samples, Axis.1, Axis.2)
write.table(ord_tab, "PCoA_TAXAbased_samples_relabund_allsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}, error=function(e){})

########################
#
# Convex-hull analysis
#
########################
setwd(paste0(as.character(args[1]), "/Figures/Ordination", sep = ""))

tryCatch({
  message("MARKER Convex hulls - phylo11 - NMDS")
  ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(Marker, NMDS1, NMDS2)
  ord_tab$NMDS1 <- round(ord_tab$NMDS1,5)
  ord_tab$NMDS2 <- round(ord_tab$NMDS2,5)
  groupNames <- as.vector(unique(ord_tab$Marker))
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
  write.table(window_area_plot, "ConvexHullAnalysis_Markers_NMDS_TAXAbased_relabund_allsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  
  message("MARKER Convex hulls - phylo11 - PCoA")
  ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, type="samples")
  ord_tab <- ord_plot$data
  ord_tab <- ord_tab %>% select(Marker, Axis.1, Axis.2)
  ord_tab$Axis.1 <- round(ord_tab$Axis.1,5)
  ord_tab$Axis.2 <- round(ord_tab$Axis.2,5)
  groupNames <- as.vector(unique(ord_tab$Marker))
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
  write.table(window_area_plot, "ConvexHullAnalysis_Markers_PCoA_TAXAbased_relabund_allsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
}, error=function(e){})


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
  write.table(window_area_plot, "ConvexHullAnalysis_replicates_NMDS_TAXAbased_relabund_allsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  
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
  write.table(window_area_plot, "ConvexHullAnalysis_replicates_PCoA_TAXAbased_relabund_allsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
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
  write.table(window_area_plot, "ConvexHullAnalysis_sites_NMDS_TAXAbased_relabund_allsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
  
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
  write.table(window_area_plot, "ConvexHullAnalysis_sites_PCoA_TAXAbased_relabund_allsamples_alltaxa.txt", sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
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
      write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_NMDS_TAXAbased_relabund_allsamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
      
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
      write.table(window_area_plot, paste0("ConvexHullAnalysis_group",i,"_PCoA_TAXAbased_relabund_allsamples_alltaxa.txt"), sep="\t", quote=F, col.names=TRUE, row.names = FALSE)
    }
    }, error=function(e){})
  }
}

########################
#
# Taxa-based ORDINATION
#
########################
#Taxa-based relative abundance (phylo11)
setwd(paste0(as.character(args[1]), "/Figures/Ordination", sep = ""))

tryCatch({
ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                            type="taxa", 
                            color = "Phylum",
                            title=paste0("NMDS - stress ",round(phylo11_stress, 3)))
ord_plot_legend <- get_legend(ord_plot)
ord_plot <- ord_plot + theme(legend.position='none')
pdf(file='NMDS_TAXAbased_taxa_relabund_allsamples_alltaxa_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='NMDS_TAXAbased_taxa_relabund_allsamples_alltaxa.pdf', width = 8, height = 8)
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
pdf(file='NMDS_TAXAbased_taxa_relabund_allsamples_alltaxa_facetPhylum_colorGenus_legend.pdf', width = 20, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='NMDS_TAXAbased_taxa_relabund_allsamples_alltaxa_facetPhylum_colorGenus.pdf', width = 15, height = 15)
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
pdf(file='PCoA_TAXAbased_taxa_relabund_allsamples_alltaxa_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='PCoA_TAXAbased_taxa_relabund_allsamples_alltaxa.pdf', width = 8, height = 8)
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
pdf(file='PCoA_TAXAbased_taxa_relabund_allsamples_alltaxa_facetPhylum_colorGenus_legend.pdf', width = 20, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='PCoA_TAXAbased_taxa_relabund_allsamples_alltaxa_facetPhylum_colorGenus.pdf', width = 15, height = 15)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                            type="samples", 
                            title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
  geom_point(size=3) +
  geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
pdf(file='NMDS_TAXAbased_samples_relabund_allsamples_alltaxa.pdf', width = 8, height = 8)
print(ord_plot)
dev.off()
}, error=function(e){})

tryCatch({
ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                            type="samples", 
                            title="PCoA") +
  geom_point(size=3) +
  geom_text_repel(aes(label = row.names(ord_plot$data)), box.padding = 0.4, min.segment.length = 0.2)
pdf(file='PCoA_TAXAbased_samples_relabund_allsamples_alltaxa.pdf', width = 8, height = 8)
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
pdf(file='NMDS_TAXAbased_taxa_relabund_allsamples_alltaxa_biplot_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='NMDS_TAXAbased_taxa_relabund_allsamples_alltaxa_biplot.pdf', width = 8, height = 8)
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
pdf(file='PCoA_TAXAbased_taxa_relabund_allsamples_alltaxa_biplot_legend.pdf', width = 11, height = 6)
print(as_ggplot(ord_plot_legend))
dev.off()
pdf(file='PCoA_TAXAbased_taxa_relabund_allsamples_alltaxa_biplot.pdf', width = 8, height = 8)
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
      pdf(file=paste0('NMDS_TAXAbased_samples_relabund_allsamples_alltaxa_group',i,'Colored_encircleReplicates.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('PCoA_TAXAbased_samples_relabund_allsamples_alltaxa_group',i,'Colored_encircleReplicates.pdf'), width = 8, height = 8)
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
      pdf(file=paste0('NMDS_TAXAbased_samples_relabund_allsamples_alltaxa_group',i,'Colored_encircleSites.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = sites, group=sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
      pdf(file=paste0('PCoA_TAXAbased_samples_relabund_allsamples_alltaxa_group',i,'Colored_encircleSites.pdf'), width = 8, height = 8)
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
      pdf(file=paste0('NMDS_TAXAbased_samples_relabund_allsamples_alltaxa_group',i,'Colored.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3)
      pdf(file=paste0('PCoA_TAXAbased_samples_relabund_allsamples_alltaxa_group',i,'Colored.pdf'), width = 8, height = 8)
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
      pdf(file=paste0('NMDS_TAXAbased_samples_relabund_allsamples_alltaxa_group',i,'Colored_encircleGroup.pdf'), width = 8, height = 8)
      print(ord_plot)
      dev.off()
      
      ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                                  type="samples", 
                                  color=paste0("group",i),
                                  title="PCoA") +
        geom_point(size=3) +
        geom_encircle(aes(fill = eval(parse(text = paste0("group",i))), group=eval(parse(text = paste0("group",i)))), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001) +
        guides(fill = FALSE)
      pdf(file=paste0('PCoA_TAXAbased_samples_relabund_allsamples_alltaxa_group',i,'Colored_encircleGroup.pdf'), width = 8, height = 8)
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
  pdf(file='NMDS_TAXAbased_samples_relabund_allsamples_alltaxa_encircleReplicates.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  
  ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                              type="samples", 
                              color="replicates",
                              title="PCoA") +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='PCoA_TAXAbased_samples_relabund_allsamples_alltaxa_encircleReplicates.pdf', width = 8, height = 8)
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
  pdf(file='NMDS_TAXAbased_samples_relabund_allsamples_alltaxa_encircleSites.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  
  ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                              type="samples", 
                              color="sites",
                              title="PCoA") +
    geom_point(size=3) +
    geom_encircle(aes(fill = sites, group = sites), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='PCoA_TAXAbased_samples_relabund_allsamples_alltaxa_encircleSites.pdf', width = 8, height = 8)
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
  pdf(file='NMDS_TAXAbased_samples_relabund_allsamples_alltaxa_encircleshapeReplicates_colorSites.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
  
  ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                              type="samples", 
                              color="sites",
                              shape="replicates",
                              title="PCoA") +
    geom_point(size=3) +
    geom_encircle(aes(fill = replicates, group=replicates), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
  pdf(file='PCoA_TAXAbased_samples_relabund_allsamples_alltaxa_encircleshapeReplicates_colorSites.pdf', width = 8, height = 8)
  print(ord_plot)
  dev.off()
}
}, error=function(e){})











tryCatch({
    ord_plot <- plot_ordination(phylo11, phylo11_NMDS.ord, 
                                type="samples",
                                color="Marker",
                                title=paste0("NMDS - stress ",round(phylo11_stress, 3))) +
      geom_point(size=3) +
      geom_encircle(aes(fill = Marker, group=Marker), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
    pdf(file='NMDS_TAXAbased_samples_relabund_allsamples_alltaxa_encircleMarkers.pdf', width = 8, height = 8)
    print(ord_plot)
    dev.off()
    
    ord_plot <- plot_ordination(phylo11, phylo11_PCoA.ord, 
                                type="samples", 
                                color="Marker",
                                title="PCoA") +
      geom_point(size=3) +
      geom_encircle(aes(fill = Marker, group=Marker), inherit.aes = TRUE, alpha = 0.4, expand = 0, s_shape = 1, spread = 0.001)
    pdf(file='PCoA_TAXAbased_samples_relabund_allsamples_alltaxa_encircleMarkers.pdf', width = 8, height = 8)
    print(ord_plot)
    dev.off()
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

setwd(paste0(as.character(args[1]), "/Figures/Network", sep = ""))

for (i in 1:9) {
  tryCatch({
    taxa_net <- make_network(filter_taxa(phylo11, function(x) sum(x >= 1) > 0, TRUE), type="taxa", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
    net_plot <- plot_network(taxa_net, filter_taxa(phylo11, function(x) sum(x >= 1) > 0, TRUE))
    pdf(file=paste0('Network_TAXAbased_taxa_relabund_allsamples_gte1perctaxa_dist0.',i,'_ASVLabeled.pdf'), width = 8, height = 8)
    print(net_plot)
    dev.off()
  }, error=function(e){})
}

for (i in 1:9) {
  tryCatch({
    sample_net <- make_network(phylo11, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
    net_plot <- plot_network(sample_net, phylo11)
    pdf(file=paste0('Network_TAXAbased_samples_relabund_allsamples_gte1perctaxa_dist0.',i,'_SampleLabeled.pdf'), width = 8, height = 8)
    print(net_plot)
    dev.off()
  }, error=function(e){})
}

if (replicateFlag == TRUE) {
  for (i in 1:9) {
    tryCatch({
      sample_net <- make_network(phylo11, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
      net_plot <- plot_network(sample_net, phylo11, color="replicates")
      pdf(file=paste0('Network_TAXAbased_samples_relabund_allsamples_gte1perctaxa_dist0.',i,'_SampleLabeled_colorReplicates.pdf'), width = 8, height = 8)
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
      pdf(file=paste0('Network_TAXAbased_samples_relabund_allsamples_gte1perctaxa_dist0.',i,'_SampleLabeled_colorSites.pdf'), width = 8, height = 8)
      print(net_plot)
      dev.off()
    }, error=function(e){})
  }
}

for (i in 1:9) {
  tryCatch({
    sample_net <- make_network(phylo11, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
    net_plot <- plot_network(sample_net, phylo11, color="Marker")
    pdf(file=paste0('Network_TAXAbased_samples_relabund_allsamples_gte1perctaxa_dist0.',i,'_SampleLabeled_colorMarkers.pdf'), width = 8, height = 8)
    print(net_plot)
    dev.off()
  }, error=function(e){})
}

########################################################################################################################
#
#
# WARNING: EVERYTHING PAST HERE REQUIRES MODIFICATIONS OF ORIGINAL phylo OBJECTS
#
#
########################################################################################################################

########################
#
# USER SPECIFIED TAXA OF INTEREST
#
########################

##################################
#
#  SUBSET phylo objects to user interest
#
##################################
if (taxaofinterestFlag == TRUE) {
  taxa_of_interest <- read.delim(as.character(args[11]), header=FALSE, stringsAsFactors=FALSE)
  taxa_of_interest <- as.vector(taxa_of_interest)
  hierarchy_level <- as.character(args[10])
  
  phylo11 <- subset_taxa(phylo11, eval(parse(text = hierarchy_level)) %in% taxa_of_interest$V1)
  
  ########################
  #
  # NETWORK ANALYSIS
  #
  ########################
  
  setwd(paste0(as.character(args[1]), "/Figures/Network", sep = ""))
  
  for (i in 1:9) {
    tryCatch({
      taxa_net <- make_network(filter_taxa(phylo11, function(x) sum(x >= 1) > 0, TRUE), type="taxa", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
      net_plot <- plot_network(taxa_net, filter_taxa(phylo11, function(x) sum(x >= 1) > 0, TRUE))
      pdf(file=paste0('Network_TAXAbased_taxa_relabund_allsamples_UserTaxaOfInterest_gte1perctaxa_dist0.',i,'_ASVLabeled.pdf'), width = 8, height = 8)
      print(net_plot)
      dev.off()
    }, error=function(e){})
  }
  
  for (i in 1:9) {
    tryCatch({
      sample_net <- make_network(phylo11, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
      net_plot <- plot_network(sample_net, phylo11)
      pdf(file=paste0('Network_TAXAbased_samples_relabund_allsamples_UserTaxaOfInterest_gte1perctaxa_dist0.',i,'_SampleLabeled.pdf'), width = 8, height = 8)
      print(net_plot)
      dev.off()
    }, error=function(e){})
  }
  
  if (replicateFlag == TRUE) {
    for (i in 1:9) {
      tryCatch({
        sample_net <- make_network(phylo11, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
        net_plot <- plot_network(sample_net, phylo11, color="replicates")
        pdf(file=paste0('Network_TAXAbased_samples_relabund_filtsamples_UserTaxaOfInterest_gte1perctaxa_dist0.',i,'_SampleLabeled_colorReplicates.pdf'), width = 8, height = 8)
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
        pdf(file=paste0('Network_TAXAbased_samples_relabund_filtsamples_UserTaxaOfInterest_gte1perctaxa_dist0.',i,'_SampleLabeled_colorSites.pdf'), width = 8, height = 8)
        print(net_plot)
        dev.off()
      }, error=function(e){})
    }
  }
  
    for (i in 1:9) {
      tryCatch({
        sample_net <- make_network(phylo11, type="samples", distance="bray", max.dist=eval(parse(text = paste0("0.",i))))
        net_plot <- plot_network(sample_net, phylo11, color="Marker")
        pdf(file=paste0('Network_TAXAbased_samples_relabund_filtsamples_UserTaxaOfInterest_gte1perctaxa_dist0.',i,'_SampleLabeled_colorMarkers.pdf'), width = 8, height = 8)
        print(net_plot)
        dev.off()
      }, error=function(e){})
    }
}
