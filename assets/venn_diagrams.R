#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

########################################
#TEMP WHILE WORKING ON SCRIPT
#args[1]<-"/Users/mcallister/Desktop/Illumina/1_WOAC/WOAC/metapipe_run/noENVblast/marker_comparison_RESULTS/NEW_script_test/comp_out" #COMP OUT directory
#args[2]<-"3" #Count number of markers
#args[3]<-"Species" #The taxonomic level to be considered (Kingdom, Phylum, Class, Order, Family, Genus, Species)
#args[4]<-"Total" #"Total" or "Terminal" Unique Taxa
########################################
library("VennDiagram")
library("dplyr")

setwd(as.character(args[1]))

numMarkers <- as.numeric(args[2])
hierarchyLevel <- as.character(args[3])
uniqueTotalTerminal <- as.character(args[4])

ind_marker <- read.delim(paste0(as.character(args[1]),"/Tables/taxacounts_individualMarkersOnly.txt", sep = ""), header=TRUE, stringsAsFactors=FALSE)
row.names(ind_marker) <- ind_marker$TaxaLevel_counts
ind_marker <- ind_marker %>% select(-TaxaLevel_counts)
keep_row <- c(paste0("unique",args[4],"_at_",args[3], sep = ""))
ind_marker <- ind_marker[rownames(ind_marker) %in% keep_row, ]

if (uniqueTotalTerminal == "Total") {
  markercompare <- read.delim(paste0(as.character(args[1]),"/Tables/taxacounts_comparisonsUniqueTotal.txt", sep = ""), header=TRUE, stringsAsFactors=FALSE)
}
if (uniqueTotalTerminal == "Terminal") {
  markercompare <- read.delim(paste0(as.character(args[1]),"/Tables/taxacounts_comparisonsUniqueTerminal.txt", sep = ""), header=TRUE, stringsAsFactors=FALSE)
}

row.names(markercompare) <- markercompare[ ,1]
markercompare <- markercompare %>% select(-starts_with("Comparison_NumUnique"))
markercompare <- markercompare %>% select(paste0(hierarchyLevel))

##################################
#
#  Make Venn Diagrams for comparisons of 2-5 markers
#
##################################
setwd(paste0(as.character(args[1]), "/Figures/VennDiagrams", sep = ""))

if (numMarkers == 2) {
  marker1 <- as.character(colnames(ind_marker)[1])
  marker2 <- as.character(colnames(ind_marker)[2])
  area1 <- as.numeric(ind_marker[1,1])
  area2 <- as.numeric(ind_marker[1,2])
  
  n12 <- as.numeric(markercompare[paste0(marker1,";",marker2, sep = ""),1])
  
  venn_plot <- draw.pairwise.venn(
    area1 = area1, area2 = area2,
    cross.area = n12,
    category = c(paste0(marker1), paste0(marker2)),
    fill = c("dodgerblue", "goldenrod1"),
    cat.col = c("dodgerblue", "goldenrod1"),
    scaled = FALSE,
    ind = TRUE
  )
  
  pdf(file=paste0('vennDiagram_2markers_',hierarchyLevel,'_unique',uniqueTotalTerminal,'.pdf', sep = ""), width = 8, height = 8)
  grid.draw(venn_plot)
  dev.off()
}

if (numMarkers == 3) {
  marker1 <- as.character(colnames(ind_marker)[1])
  marker2 <- as.character(colnames(ind_marker)[2])
  marker3 <- as.character(colnames(ind_marker)[3])
  area1 <- as.numeric(ind_marker[1,1])
  area2 <- as.numeric(ind_marker[1,2])
  area3 <- as.numeric(ind_marker[1,3])
  
  num12only <- as.numeric(markercompare[paste0(marker1,";",marker2, sep = ""),1])
  num13only <- as.numeric(markercompare[paste0(marker1,";",marker3, sep = ""),1])
  num23only <- as.numeric(markercompare[paste0(marker2,";",marker3, sep = ""),1])
  num123only <- as.numeric(markercompare[paste0(marker1,";",marker2,";",marker3, sep = ""),1])
  
  n12 <- num12only + num123only
  n13 <- num13only + num123only
  n23 <- num23only + num123only
  n123 <- num123only
  
  venn_plot <- draw.triple.venn(
    area1 = area1, area2 = area2, area3 = area3,
    n12 = n12, n13 = n13, n23 = n23, n123 = n123,
    category = c(paste0(marker1), paste0(marker2), paste0(marker3)),
    fill = c("dodgerblue", "goldenrod1", "seagreen3"),
    cat.col = c("dodgerblue", "goldenrod1", "seagreen3"),
    scaled = FALSE,
    ind = TRUE
  )
  
  pdf(file=paste0('vennDiagram_3markers_',hierarchyLevel,'_unique',uniqueTotalTerminal,'.pdf', sep = ""), width = 8, height = 8)
  grid.draw(venn_plot)
  dev.off()
}

if (numMarkers == 4) {
  marker1 <- as.character(colnames(ind_marker)[1])
  marker2 <- as.character(colnames(ind_marker)[2])
  marker3 <- as.character(colnames(ind_marker)[3])
  marker4 <- as.character(colnames(ind_marker)[4])
  area1 <- as.numeric(ind_marker[1,1])
  area2 <- as.numeric(ind_marker[1,2])
  area3 <- as.numeric(ind_marker[1,3])
  area4 <- as.numeric(ind_marker[1,4])
  
  num12only <- as.numeric(markercompare[paste0(marker1,";",marker2, sep = ""),1])
  num13only <- as.numeric(markercompare[paste0(marker1,";",marker3, sep = ""),1])
  num14only <- as.numeric(markercompare[paste0(marker1,";",marker4, sep = ""),1])
  num23only <- as.numeric(markercompare[paste0(marker2,";",marker3, sep = ""),1])
  num24only <- as.numeric(markercompare[paste0(marker2,";",marker4, sep = ""),1])
  num34only <- as.numeric(markercompare[paste0(marker3,";",marker4, sep = ""),1])
  num123only <- as.numeric(markercompare[paste0(marker1,";",marker2,";",marker3, sep = ""),1])
  num124only <- as.numeric(markercompare[paste0(marker1,";",marker2,";",marker4, sep = ""),1])
  num134only <- as.numeric(markercompare[paste0(marker1,";",marker3,";",marker4, sep = ""),1])
  num234only <- as.numeric(markercompare[paste0(marker2,";",marker3,";",marker4, sep = ""),1])
  num1234only <- as.numeric(markercompare[paste0(marker1,";",marker2,";",marker3,";",marker4, sep = ""),1])
  
  n12 <- num12only + num123only + num124only + num1234only
  n13 <- num13only + num123only + num134only + num1234only
  n14 <- num14only + num124only + num134only + num1234only
  n23 <- num23only + num123only + num234only + num1234only
  n24 <- num24only + num124only + num234only + num1234only
  n34 <- num34only + num134only + num234only + num1234only
  n123 <- num123only + num1234only
  n124 <- num124only + num1234only
  n134 <- num134only + num1234only
  n234 <- num234only + num1234only
  n1234 <- num1234only
  
  venn_plot <- draw.quad.venn(
    area1 = area1, area2 = area2, area3 = area3, area4 = area4,
    n12 = n12, n13 = n13, n14 = n14, n23 = n23, n24 = n24, n34 = n34,
    n123 = n123, n124 = n124, n134 = n134, n234 = n234,
    n1234 = n1234,
    category = c(paste0(marker1), paste0(marker2), paste0(marker3), paste0(marker4)),
    fill = c("dodgerblue", "goldenrod1", "seagreen3", "darkorange1"),
    cat.col = c("dodgerblue", "goldenrod1", "seagreen3", "darkorange1"),
    scaled = FALSE,
    ind = TRUE
  )
  
  pdf(file=paste0('vennDiagram_4markers_',hierarchyLevel,'_unique',uniqueTotalTerminal,'.pdf', sep = ""), width = 8, height = 8)
  grid.draw(venn_plot)
  dev.off()
}

if (numMarkers == 5) {
  marker1 <- as.character(colnames(ind_marker)[1])
  marker2 <- as.character(colnames(ind_marker)[2])
  marker3 <- as.character(colnames(ind_marker)[3])
  marker4 <- as.character(colnames(ind_marker)[4])
  marker5 <- as.character(colnames(ind_marker)[5])
  area1 <- as.numeric(ind_marker[1,1])
  area2 <- as.numeric(ind_marker[1,2])
  area3 <- as.numeric(ind_marker[1,3])
  area4 <- as.numeric(ind_marker[1,4])
  area5 <- as.numeric(ind_marker[1,5])
  
  num12only <- as.numeric(markercompare[paste0(marker1,";",marker2, sep = ""),1])
  num13only <- as.numeric(markercompare[paste0(marker1,";",marker3, sep = ""),1])
  num14only <- as.numeric(markercompare[paste0(marker1,";",marker4, sep = ""),1])
  num15only <- as.numeric(markercompare[paste0(marker1,";",marker5, sep = ""),1])
  num23only <- as.numeric(markercompare[paste0(marker2,";",marker3, sep = ""),1])
  num24only <- as.numeric(markercompare[paste0(marker2,";",marker4, sep = ""),1])
  num25only <- as.numeric(markercompare[paste0(marker2,";",marker5, sep = ""),1])
  num34only <- as.numeric(markercompare[paste0(marker3,";",marker4, sep = ""),1])
  num35only <- as.numeric(markercompare[paste0(marker3,";",marker5, sep = ""),1])
  num45only <- as.numeric(markercompare[paste0(marker4,";",marker5, sep = ""),1])
  num123only <- as.numeric(markercompare[paste0(marker1,";",marker2,";",marker3, sep = ""),1])
  num124only <- as.numeric(markercompare[paste0(marker1,";",marker2,";",marker4, sep = ""),1])
  num125only <- as.numeric(markercompare[paste0(marker1,";",marker2,";",marker5, sep = ""),1])
  num134only <- as.numeric(markercompare[paste0(marker1,";",marker3,";",marker4, sep = ""),1])
  num135only <- as.numeric(markercompare[paste0(marker1,";",marker3,";",marker5, sep = ""),1])
  num145only <- as.numeric(markercompare[paste0(marker1,";",marker4,";",marker5, sep = ""),1])
  num234only <- as.numeric(markercompare[paste0(marker2,";",marker3,";",marker4, sep = ""),1])
  num235only <- as.numeric(markercompare[paste0(marker2,";",marker3,";",marker5, sep = ""),1])
  num245only <- as.numeric(markercompare[paste0(marker2,";",marker4,";",marker5, sep = ""),1])
  num345only <- as.numeric(markercompare[paste0(marker3,";",marker4,";",marker5, sep = ""),1])
  num1234only <- as.numeric(markercompare[paste0(marker1,";",marker2,";",marker3,";",marker4, sep = ""),1])
  num1235only <- as.numeric(markercompare[paste0(marker1,";",marker2,";",marker3,";",marker5, sep = ""),1])
  num1245only <- as.numeric(markercompare[paste0(marker1,";",marker2,";",marker4,";",marker5, sep = ""),1])
  num1345only <- as.numeric(markercompare[paste0(marker1,";",marker3,";",marker4,";",marker5, sep = ""),1])
  num2345only <- as.numeric(markercompare[paste0(marker2,";",marker3,";",marker4,";",marker5, sep = ""),1])
  num12345only <- as.numeric(markercompare[paste0(marker1,";",marker2,";",marker3,";",marker4,";",marker5, sep = ""),1])
  
  n12 <- num12only + num123only + num124only + num125only + num1234only + num1235only + num1245only + num12345only
  n13 <- num13only + num123only + num134only + num135only + num1234only + num1235only + num1345only + num12345only
  n14 <- num14only + num124only + num134only + num145only + num1234only + num1245only + num1345only + num12345only
  n15 <- num15only + num125only + num135only + num145only + num1235only + num1245only + num1345only + num12345only
  n23 <- num23only + num123only + num234only + num235only + num1234only + num1235only + num2345only + num12345only
  n24 <- num24only + num124only + num234only + num245only + num1234only + num1245only + num2345only + num12345only
  n25 <- num25only + num125only + num235only + num245only + num1235only + num1245only + num2345only + num12345only
  n34 <- num34only + num134only + num234only + num345only + num1234only + num1345only + num2345only + num12345only
  n35 <- num35only + num135only + num235only + num345only + num1235only + num1345only + num2345only + num12345only
  n45 <- num45only + num145only + num245only + num345only + num1245only + num1345only + num2345only + num12345only
  n123 <- num123only + num1234only + num1235only + num12345only
  n124 <- num124only + num1234only + num1245only + num12345only
  n125 <- num125only + num1235only + num1245only + num12345only
  n134 <- num134only + num1234only + num1345only + num12345only
  n135 <- num135only + num1235only + num1345only + num12345only
  n145 <- num145only + num1245only + num1345only + num12345only
  n234 <- num234only + num1234only + num2345only + num12345only
  n235 <- num235only + num1235only + num2345only + num12345only
  n245 <- num245only + num1245only + num2345only + num12345only
  n345 <- num345only + num1345only + num2345only + num12345only
  n1234 <- num1234only + num12345only
  n1235 <- num1235only + num12345only
  n1245 <- num1245only + num12345only
  n1345 <- num1345only + num12345only
  n2345 <- num2345only + num12345only
  n12345 <- num12345only
  
  venn_plot <- draw.quintuple.venn(
    area1 = area1, area2 = area2, area3 = area3, area4 = area4, area5 = area5,
    n12 = n12, n13 = n13, n14 = n14, n15 = n15,
    n23 = n23, n24 = n24, n25 = n25,
    n34 = n34, n35 = n35, n45 = n45,
    n123 = n123, n124 = n124, n125 = n125, n134 = n134, n135 = n135, n145 = n145, n234 = n234, n235 = n235, n245 = n245, n345 = n345,
    n1234 = n1234, n1235 = n1235, n1245 = n1245, n1345 = n1345, n2345 = n2345, n12345 = n12345,
    category = c(paste0(marker1), paste0(marker2), paste0(marker3), paste0(marker4), paste0(marker5)),
    fill = c("dodgerblue", "goldenrod1", "seagreen3", "darkorange1", "orchid3"),
    cat.col = c("dodgerblue", "goldenrod1", "seagreen3", "darkorange1", "orchid3"),
    scaled = FALSE,
    ind = TRUE
  )
  
  pdf(file=paste0('vennDiagram_5markers_',hierarchyLevel,'_unique',uniqueTotalTerminal,'.pdf', sep = ""), width = 8, height = 8)
  grid.draw(venn_plot)
  dev.off()
}
