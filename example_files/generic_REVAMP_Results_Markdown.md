# REVAMP Results Markdown
---

#### **Created by: 'Omics Group, NOAA PMEL/UW CICOES**

## Introduction

**Text Describing what REVAMP is for.** Lorem ipsum dolor sit amet, consectetur adipiscing elit. Curabitur velit nisl, condimentum a enim quis, facilisis imperdiet augue. Vivamus suscipit aliquam felis pharetra ultricies. Interdum et malesuada fames ac ante ipsum primis in faucibus. Curabitur pharetra, augue a porta laoreet, tellus elit viverra arcu, ac vehicula felis metus nec eros. Ut tempus, nisl non malesuada venenatis, dui est viverra turpis, sed iaculis lectus lacus id ante. Mauris ac dui at risus auctor sollicitudin ac ac erat. Nulla tempor enim ac ipsum consectetur vestibulum. Phasellus vel massa eget dolor consectetur commodo a vitae neque.

Duis vulputate felis sit amet elit pretium efficitur. Ut accumsan libero sed convallis cursus. Donec posuere ultrices dolor, sollicitudin molestie ante. Nulla quis quam in risus facilisis viverra. Proin suscipit sit amet magna eu mollis. Aenean rutrum ac nunc at aliquet. Nullam velit nisi, dignissim a interdum eget, hendrerit sit amet magna. Morbi leo nisi, aliquet id velit vel, malesuada pellentesque mauris. Maecenas interdum, felis nec ornare vestibulum, orci justo vestibulum nulla, non eleifend magna velit ullamcorper orci. Praesent sed facilisis libero. Interdum et malesuada fames ac ante ipsum primis in faucibus. Aenean venenatis ipsum vel enim facilisis, vitae ornare libero vulputate. Donec eu massa pretium, auctor tortor ut, pretium mi. Vestibulum varius mi non porta faucibus. Vestibulum scelerisque diam vulputate tortor pellentesque, ac aliquet mi ullamcorper.

If you use REVAMP, please cite:

```
REVAMP citation. In preparation.
```

In addition, REVAMP uses some third-party software that should also be cited:

```
dada2 (including blog post)
cutadapt
blast
Krona
taxonkit
phyloseq
vegan

```

---

## REVAMP Directory Structure

```
OUTDIR
├── run_logs
├── cutadapt
├── dada2
├── blast_results
├── ASV2Taxonomy
├── processed_tables
│   └── replicate_based_detection (if replicates given)
└── Figures
    ├── 00_KRONA_plots
    ├── 01_Maps
    ├── 02_Barcharts
    │   └── read_count & relative_abundance
    ├── 03_Heatmaps
    │   └── ASV_based & Taxonomy_merge_based
    ├── 04_Alpha_Diversity
    │   └── ASV_based & Taxonomy_merge_based
    ├── 05_Ordination
    │   ├── ASV_based: read_count & relative_abundance
    │   └── Taxonomy_merge_based
    │       ├── filterInclude_TOSPECIES_only
    │       │   └── read_count & relative_abundance
    │       └── read_count & relative_abundance
    ├── 06_Network
    │   ├── ASV_based: read_count & relative_abundance
    │   └── Taxonomy_merge_based: read_count & relative_abundance
    ├── 07_Rarefaction_Curves
    ├── 08_EnvironmentFit_Ordination
    │   └── ASV_based & Taxonomy_merge_based
    ├── ReadsVSReplicateDetection (if replicates given)
    └── Taxa_of_interest (if given)
        ├── 02_Barcharts
        ├── 03_Heatmaps
        └── 06_Network

```

All information printed to screen during a run is stored to `run.log`. If more than one run is done on the same folder, the directory `run_logs` is created and the previous `run.log` is renamed to append the date/time and moved to that directory.

Raw reads are first filtered using `cutadapt`, to identify and trim the PCR primers. This does enable partitioning of two PCR products from the same barcode. Then, trimmed reads are brought into `dada2`, where they are quality trimmed and filtered (step 1) followed by learning error, dereplication, merging of F/R reads, and determination of ASVs (Amplicon Sequence Variants) (step 2). ASVs are then "blasted" using BLASTn (`blast_results` directory) against NCBI's *nt*. These results are parsed, assigned to taxonomy (`ASV2Taxonomy` directory), and various table products (`processed_tables`) and `Figures` are generated.

The `Figure` directory includes four primary choices for preference of data analysis (where applicable). The User can choose to look at figures where ASVs are considered a single ecological unit (`ASV_based`) or where ASVs have been collapsed so that each represents a single unique taxonomic hierarchy (`Taxonomy_merge_based`). Further, figures can be generated based on the raw count data (`read_count`) or based on the normalized relative percent abundance (`relative_abundance`).

## REVAMP Files of Interest

All configuration files are copied to the REVAMP results folder. `config_file.txt` controls all of the main options at the front end of REVAMP. `figure_config_file.txt` controls all the options for the back end figure generation of REVAMP. The provided sample metadata file is copied (`sample_metadata.txt`) and converted/cleaned for use in R (`sample_metadata_forR.txt`).

If you wish to start and stop REVAMP at varying points along the pipeline, this can be done. The `progress.txt` file is used to track the completion of each checkpoint. If you wish to redo a checkpoint that has already passed, simply delete sequentially from the bottom of the `progress.txt` file until you have deleted the step you want to redo. **It is not recommended to delete a step in the middle while allowing the other steps to remain. It will break.**

Each R script outputs stdout and stderror to a log file within the R scripts primary directory. In addition, if you wish to customize or debug any R-based figure/file, you can open any of the R scripts in the `assets` folder. The second field of the tab-delimited file `Rscript_arguments.log` gives each of the positional arguments fed to each script during the pipeline. Simply uncomment the args block at the front end of the script and input each argument.

Besides the obvious outputs in the `processed_tables` and `Figures` directory, the User may also find these files useful:

* `OUTDIR/dada2/ASVs.fa` - These are the fasta sequences used in all subsequent steps.
* `OUTDIR/dada2/ASVs_counts.tsv` - This biom table is the original table showing read counts per ASV per sample.
* `OUTDIR/blast_results/ASV_blastn_nt_formatted.txt` - Shows the list of identical best hits for each OTU, with the taxids listed.
* `OUTDIR/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy.tsv` - This biom table shows the read counts per ASV per sample, where unique taxonomy is collapsed (summed) to be represented by only a single ASV.
* The ASVs_counts... files are also represented in `OUTDIR/ASV2Taxonomy` with "Unknowns" removed.
* `OUTDIR/ASV2Taxonomy/OUTDIR_unknown_asvids.txt` can help to track the reasons for "Unknown" assignments.
* `OUTDIR/ASV2Taxonomy/OUTDIR_heatmap_multiASV.txt` Shows a heatmap of read counts for taxonomic assignments with more than one ASV. This kind of result could show where multi-ASVs with the same taxonomy might be the result of erroneous sequence/chimeric error (i.e. one ASV with lots of reads with a few others with the same taxonomy but very few sporadic read hits) as well as which multi-ASVs might represent different ecotypes of the same taxa (i.e. multiple ASVs with believable read distrubution among samples).


## User configuration

```
Output directory: $outdirectory
Input Read Directory: $readfolderpath

Cutadapt parameters:
   Forward Primer: $primerF
   Reverse Primer: $primerR
Dada2 parameters:
   Failed merge? Use: $failedMerge_useDirection
   Min length: $dada_minlength
   Remove PhiX? $dada_phix
   TrunQ: $dada_trunQ
   MaxEE1: $dada_maxEE1
   MaxEE2: $dada_maxEE2
   Trim from Right: $dada_trimRight
   Trim from Left: $dada_trimLeft
BLAST parameters:
   Location of nt: $locationNTdatabase
   BLAST length cutoff: $blastLengthCutoff
   BLAST mode: $blastMode
Assign taxonomy parameters:
   Taxonomy confidence cutoffs (P,C,O,F,G,S): $taxonomyCutoffs
   File designating ASVs to remove: $removeASVsFILE
Interpreted sample parameters:
   Controls Present: $controlspresent
   Positive Controls: $controlPos
   Negative Controls: $controlNeg
   Replicates: $replicates
   Sites: $sites
   Groups present: $groupsDefinedFlag
   Number of groups: $numberGroupsDefined
   Chemistry columns: $chemData
Figure parameters:
   Filter low quality samples: $filterLowQualSamples
   Filter threshold for low qual: $filterPercentLowQualSamples
   Filter threshold %, < to zzOther: $filterPercent
   Remove NAs from figures: $removeNA
   Taxa of interest? $providedTaxaOfInterest
   Taxa of interest hierarchy level: $taxaOfInterestLevel
   Scale of pie charts on map: $pieScale
```

## Result Stats

```
Number of Samples: $numSamples
Number of ASVs: $numASVs
Number of "ASVs" when collapsed on taxonomy: $numASVsTaxonomyMerge
Average Length: $avgLength
Average No. Reads In: $avgNumReadsIN
Average No. Reads Out: $avgNumReadsOUT
Average Percent Reads Passing: $avgPercPass
Average Percent ASVs w/ Taxonomy Classification: $avgPercTaxClass
Total Hits Terminal at Species: $totalSpeciesHits
Total Hits Terminal at Genus: $totalGenusHits
Total Hits Terminal at Family: $totalFamilyHits
Total Hits Terminal at Order: $totalOrderHits
Total Hits Terminal at Class: $totalClassHits
Total Hits Terminal at Phylum: $totalPhylumHits
```

REPLACEsampleSTATS

## Result Figures (sampling)

### Maps

Several maps are created using the lat/long from the sample metadata file. If "sites" are given, only one data point for each site is displayed.

<figure>
<img src="/Users/mcallister/Desktop/chris_test/18S/CP_all_RUN20220208_out/Figures/01_Maps/mapBathy_sanslegend_datapoints.png" alt="bathymap" width="600"/>
<figcaption align = "center">Figure 1 - Bathymetric map of sampling area. <code>OUTDIR/Figures/01_Maps/mapBathy_sanslegend_datapoints.pdf</code> </figcaption>
</figure>

### Krona Plots

The best place to start looking at results is `OUTDIR/Figures/00_KRONA_plots`. A Krona plot is an interactive hierarchical data viewer loaded in an internet browser of your choice (Chrome recommended). You can search for taxa of interest, and zoom to selected taxa, seeing total number of reads in the top right corner. An excellent method for looking at the data quickly and effectively. Krona plots are provided by sample (`OUTDIR_master_krona.html`) or with all samples summed to see all the taxa in the reads that were processed (`OUTDIR_samplesSummedKRONA.html`).

### Bar Charts

Next, it is a good idea to look at the overall evenness of your sequencing effort (i.e. compare the number of reads passing QC between samples).

<figure>
<img src="/Users/mcallister/Desktop/chris_test/18S/CP_all_RUN20220208_out/Figures/02_Barcharts/read_count/barplot_readcount_allsamples_alltaxa_Phylum.png" alt="barchart_readcount"/>
<figcaption align = "center">Figure 2 - Bar chart showing raw read counts per sample (legend available in folder with sameName_legend). <code>OUTDIR/Figures/02_Barcharts/read_count/barplot_readcount_allsamples_alltaxa_Phylum.pdf</code> </figcaption>
</figure>

There are many options for viewing the bar charts. One easier-to-visualize method is to filter all ASVs <x% to the zzOther category. This allows us to view and interpret more highly abundant organisms. Many of these bar charts exist averaging the ASV relative abundance of samples coming from the same replicates/sites together ([replicate/site]GroupedSamples).

<figure>
<img src="/Users/mcallister/Desktop/chris_test/18S/CP_all_RUN20220208_out/Figures/02_Barcharts/relative_abundance/barplot_relabund_filtsamples_filtLT5PERCtaxa_Family.png" alt="barchart_relabund"/>
<figcaption align = "center">Figure 3A - Bar chart showing relative abundance (%) of taxa per sample. <code>OUTDIR/Figures/02_Barcharts/relative_abundance/barplot_relabund_filtsamples_filtLT5PERCtaxa_Family.pdf</code> </figcaption>
</figure>

<figure>
<img src="/Users/mcallister/Desktop/chris_test/18S/CP_all_RUN20220208_out/Figures/02_Barcharts/relative_abundance/barplot_relabund_filtsamples_filtLT5PERCtaxa_Family_legend.png" alt="barchart_relabund_legend" width="800"/>
<figcaption align = "center">Figure 3B - Bar chart legend. Legends are stripped from figures owing to their fluctuating size. <code>OUTDIR/Figures/02_Barcharts/relative_abundance/barplot_relabund_filtsamples_filtLT5PERCtaxa_Family_legend.pdf</code> </figcaption>
</figure>

### Rarefaction Curves

Rarefaction curves show the number of species found at a particular effort for each sample. This provides an indication if enough sampling has been done to capture a majority of species in the sample. It also gives a visual representation of the differences of richness between samples.

<figure>
<img src="/Users/mcallister/Desktop/chris_test/18S/CP_all_RUN20220208_out/Figures/07_Rarefaction_Curves/RarefactionCurve_ASVbased_readcount_allsamples_alltaxa.png" alt="rarefaction"/>
<figcaption align = "center">Figure 4 - Rarefaction curve showing all samples. <code>OUTDIR/Figures/07_Rarefaction_Curves/RarefactionCurve_ASVbased_readcount_allsamples_alltaxa.pdf</code> </figcaption>
</figure>

### Heatmaps

Heatmaps are provided and can reveal patterns in sample clustering.

### Alpha Deversity Metrics

Alpha diversity metrics are presented in tabular form, as well as visually - including coloring by different chemistry and groups (if provided) and sorted on Chao1 and Shannon metrics.

<figure>
<img src="/Users/mcallister/Desktop/chris_test/18S/CP_all_RUN20220208_out/Figures/04_Alpha_Diversity/Taxonomy_merge_based/alphaDiversity_Taxa_normalized_filtsamples_alltaxa_allMeasures.png" alt="alphadiv"/>
<figcaption align = "center">Figure 5 - Alpha diversity metrics, facetted by metric (Observed, Chao1, ACE, Shannon, Simpson, Inverse Simpson, and Fisher. <code>OUTDIR/Figures/04_Alpha_Diversity/Taxonomy_merge_based/alphaDiversity_Taxa_normalized_filtsamples_alltaxa_allMeasures.pdf</code> </figcaption>
</figure>

### Ordination

Ordination can be a very fast method for comparing the similarities in biodiversity between samples. REVAMP provides several ordinations, through NMDS and PCoA, colored and encircled based on the groupings, chemistry data, and replicates/sites given in the sample metadata. Note that the `Taxonomy_merge_based` folder contains a separate analysis using only those taxa assigned to the species level (`filterInclude_TOSPECIES_only`).

<figure>
<img src="/Users/mcallister/Desktop/chris_test/18S/CP_all_RUN20220208_out/Figures/05_Ordination/ASV_based/relative_abundance/NMDS_ASVbased_samples_relabund_filtsamples_alltaxa.png" alt="NMDS_samples"/>
<figcaption align = "center">Figure 6 - NMDS ordination plot showing all samples with labels. Stress value indicates goodness of fit with values above 0.2 suspect. <code>OUTDIR/Figures/05_Ordination/ASV_based/relative_abundance/NMDS_ASVbased_samples_relabund_filtsamples_alltaxa.pdf</code> </figcaption>
</figure>

IFREPLICATESORSITES@@>> (what to do if both are given)

<figure>
<img src="/Users/mcallister/Desktop/chris_test/18S/CP_all_RUN20220208_out/Figures/05_Ordination/ASV_based/relative_abundance/NMDS_ASVbased_samples_relabund_filtsamples_alltaxa_encircleReplicates.png" alt="NMDS_encircleReplicates"/>
<figcaption align = "center">Figure Replicate 1 - NMDS ordination plot showing all samples encircled by their replicate assignments. Stress value indicates goodness of fit with values above 0.2 suspect. <code>OUTDIR/Figures/05_Ordination/ASV_based/relative_abundance/NMDS_ASVbased_samples_relabund_filtsamples_alltaxa_encircleReplicates.pdf</code> </figcaption>
</figure>

<<@@

The biplot (combination of samples and taxa/ASVs) and the facetted plot (by Phylum, with Genera colored) are both useful plots for examining the relative contribution of taxa on the distribution of samples within the ordination space.

<figure>
<img src="/Users/mcallister/Desktop/chris_test/18S/CP_all_RUN20220208_out/Figures/05_Ordination/ASV_based/relative_abundance/NMDS_ASVbased_taxa_relabund_filtsamples_alltaxa_biplot.png" alt="NMDS_biplot"/>
<figcaption align = "center">Figure 7A - NMDS ordination biplot showing all samples (medium salmon dots) and ASVs (small dots). <code>OUTDIR/Figures/05_Ordination/ASV_based/relative_abundance/NMDS_ASVbased_taxa_relabund_filtsamples_alltaxa_biplot.pdf</code> </figcaption>
</figure>

<figure>
<img src="/Users/mcallister/Desktop/chris_test/18S/CP_all_RUN20220208_out/Figures/05_Ordination/ASV_based/relative_abundance/NMDS_ASVbased_taxa_relabund_filtsamples_alltaxa_biplot_legend.png" alt="NMDS_biplotleg"/>
<figcaption align = "center">Figure 7B - NMDS ordination biplot legend. <code>OUTDIR/Figures/05_Ordination/ASV_based/relative_abundance/NMDS_ASVbased_taxa_relabund_filtsamples_alltaxa_biplot_legend.pdf</code> </figcaption>
</figure>

<figure>
<img src="/Users/mcallister/Desktop/chris_test/18S/CP_all_RUN20220208_out/Figures/05_Ordination/ASV_based/relative_abundance/NMDS_ASVbased_taxa_relabund_filtsamples_alltaxa_facetPhylum_colorGenus.png" alt="NMDS_facet"/>
<figcaption align = "center">Figure 8A - NMDS ordination with Phyla facetted and colored internally by Genera. <code>OUTDIR/Figures/05_Ordination/ASV_based/relative_abundance/NMDS_ASVbased_taxa_relabund_filtsamples_alltaxa_facetPhylum_colorGenus.pdf</code> </figcaption>
</figure>

<figure>
<img src="/Users/mcallister/Desktop/chris_test/18S/CP_all_RUN20220208_out/Figures/05_Ordination/ASV_based/relative_abundance/NMDS_ASVbased_taxa_relabund_filtsamples_alltaxa_facetPhylum_colorGenus_legend.png" alt="NMDS_facetleg"/>
<figcaption align = "center">Figure 8B - NMDS ordination with Phyla facetted and colored internally by Genera. <code>OUTDIR/Figures/05_Ordination/ASV_based/relative_abundance/NMDS_ASVbased_taxa_relabund_filtsamples_alltaxa_facetPhylum_colorGenus_legend.pdf</code> </figcaption>
</figure>

### Network Analysis

Network analysis can be a valuable tool for visualizing the distance/similarities between samples and ASVs. In the folders of this analysis, networks with edges (lines) connecting nodes (samples or ASVs) are drawn with different Bray distance metric thresholds. Samples connecting at smaller thresholds are more similar. If replicates/sites are given, a colored version will be created.

<figure>
<img src="/Users/mcallister/Desktop/chris_test/18S/CP_all_RUN20220208_out/Figures/06_Network/Taxonomy_merge_based/relative_abundance/Network_TAXAbased_samples_relabund_filtsamples_gte1perctaxa_dist0.4_SampleLabeled.png" alt="network"/>
<figcaption align = "center">Figure 9 - Network diagram showing the connections between samples. Any two samples with a line connecting them meet the Bray distance metric threshold. <code>OUTDIR/Figures/06_Network/Taxonomy_merge_based/relative_abundance/Network_TAXAbased_samples_relabund_filtsamples_gte1perctaxa_dist0.4_SampleLabeled.pdf</code> </figcaption>
</figure>

### Environmental Fit Ordination

Using the `vegan` R package, REVAMP creates a new ordination plot (NMDS only) and runs environmental fitting on the continuos parameters given in the "chem" columns (i.e. non-controlled vocabulary, if given) of the sample metadata file, as well as on all the relative abundance percentages for each ASV. Output includes tables with the ordination coordinates (axes 1/2), r-squared, and p-value, as well as figures labelling all the samples in the new ordination space, encircling replicates/sites (if given), and ordination plots showing the "chem" values only (if given) and all vectors (including ASVs) at 0.001 and 0.05 p-value significance.

<figure>
<img src="/Users/mcallister/Desktop/chris_test/18S/CP_all_RUN20220208_out/Figures/08_EnvironmentFit_Ordination/Taxonomy_merge_based/NMDS_vegan_TAXAbased_samples_relabund_filtsamples_alltaxa.png" alt="vegan_NMDS"/>
<figcaption align = "center">Figure 10 - NMDS ordination plot calculated in vegan (R package) showing all samples with labels. <code>OUTDIR/Figures/08_EnvironmentFit_Ordination/Taxonomy_merge_based/NMDS_vegan_TAXAbased_samples_relabund_filtsamples_alltaxa.pdf</code> </figcaption>
</figure>

<figure>
<img src="/Users/mcallister/Desktop/chris_test/18S/CP_all_RUN20220208_out/Figures/08_EnvironmentFit_Ordination/Taxonomy_merge_based/NMDS_vegan_TAXAbased_samples_relabund_filtsamples_alltaxa_environmentFitVectors_p0.001.png" alt="vegan_NMDS_envparam"/>
<figcaption align = "center">Figure 11 - NMDS ordination plot calculated in vegan (R package), with all 0.001 p-value significant vectors shown, including chemistry (if given) and ASVs. <code>OUTDIR/Figures/08_EnvironmentFit_Ordination/Taxonomy_merge_based/NMDS_vegan_TAXAbased_samples_relabund_filtsamples_alltaxa_environmentFitVectors_p0.001.pdf</code> </figcaption>
</figure>

### Taxa-of-Interest Folder

If a file giving taxa of interest (at a given taxonomic hierarchy) is given, REVAMP will create bar charts, heatmaps, and network figures showing only ASVs matching those taxonomies.

### Replicate flag only

If a `replicates` column is provided in the sample metadata file, REVAMP will calculate some useful tables and figures. In the `processed_tables` directory, a directory named `replicate_based_detection` will be created. There, you can find tables of stats on the relative abundance of each ASV and how many replicates it was found in for each replicate group, as well as presence/absence tables suitable for importing into unmarked (R package) for performing occupancy modelling. In the `Figures` directory, `ReadsVSReplicateDetection` contains violin boxplots showing how variable ASV abundance is depending on how many replicates they are found in.

---

##### Legal Disclaimer

*This repository is a software product and is not official communication
of the National Oceanic and Atmospheric Administration (NOAA), or the
United States Department of Commerce (DOC).  All NOAA GitHub project
code is provided on an 'as is' basis and the user assumes responsibility
for its use.  Any claims against the DOC or DOC bureaus stemming from
the use of this GitHub project will be governed by all applicable Federal
law.  Any reference to specific commercial products, processes, or services
by service mark, trademark, manufacturer, or otherwise, does not constitute
or imply their endorsement, recommendation, or favoring by the DOC.
The DOC seal and logo, or the seal and logo of a DOC bureau, shall not
be used in any manner to imply endorsement of any commercial product
or activity by the DOC or the United States Government.*