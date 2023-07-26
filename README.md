![REVAMP_logo_v3_color-tagline](https://github.com/McAllister-NOAA/REVAMP/assets/60410177/00f0f999-043d-43c9-8bf8-1b8aa86bebca)

# REVAMP: Rapid Exploration and Visualization through an Automated Metabarcoding Pipeline

Created by Sean M. McAllister, Christopher Paight, Emily L. Norton, and Matthew P. Galaska

**REVAMP** is designed to streamline the processing of metabarcoding data from raw fastq data files to the generation of information and visuals. The purpose is to standardize methods for rapid assessments of molecular ecology datasets, including those for environmental DNA (eDNA) monitoring.

Developed by the **Ocean Molecular Ecology** ([OME](https://www.pmel.noaa.gov/ocean-molecular-ecology/)) Group at the Pacific Marine Environmental Laboratory ([PMEL/NOAA](https://www.pmel.noaa.gov/)) in cooperation with the Cooperative Institute for Climate, Ocean, & Ecosystem Studies ([CICOES/UW](https://cicoes.uw.edu/)).

## Citation

**If you find REVAMP useful in your research, please cite...**

```
IN REVIEW
```

**REVAMP is a wrapper for several tools, which should also be cited:**

* Cutadapt – [Martin, 2011](https://doi.org/10.14806/ej.17.1.200)
* DADA2 – [Callahan et al., 2016](https://doi.org/10.1038/nmeth.3869) ([useful tutorial](https://benjjneb.github.io/dada2/index.html))
* BLASTn – [Camacho et al., 2009](https://doi.org/10.1186/1471-2105-10-421)
* TaxonKit – [Shen and Ren, 2021](https://doi.org/10.1016/j.jgg.2021.03.006)
* Krona – [Ondov et al., 2011](https://doi.org/10.1186/1471-2105-12-385)
* phyloseq – [McMurdie and Holmes, 2013](https://doi.org/10.1371/journal.pone.0061217) ([useful tutorial](https://joey711.github.io/phyloseq/))
* vegan – [Oksanen et al., 2020](https://cran.r-project.org/web/packages/vegan/)

Many other auxiliary R packages make the figures possible: see [dependencies](https://github.com/McAllister-NOAA/REVAMP#install-the-necessary-dependencies-for-revamp).

## REVAMP Workflow


|![FlowChart Figure_web](https://github.com/McAllister-NOAA/REVAMP/assets/60410177/e4cc4c01-3f62-40bf-a951-fcfebe557bf4)|
|:--:|
|REVAMP primary workflow from inputs to table and figure products. Independent checkpoints in the pipeline are indicated with the dashed red line. Single asterisk designates where user-supplied BLASTn results are input in the pipeline. Doubled asterisk designates where alternative SILVAngs-based taxonomy assignments are supplemented for taxonomy decisions made by the pipeline.|

In the primary REVAMP workflow (above), raw demultiplexed reads from the sequencer (fastq/fastq.gz) are adaptor- and primer-trimmed using Cutadapt (v. 3.7) ([Martin, 2011](https://doi.org/10.14806/ej.17.1.200)), requiring forward or reverse primer matches for each respective read and trimming the opposing primer when detected. Next, adaptor- and primer-trimmed reads are passed to DADA2 (v 1.18) ([Callahan et al., 2016](https://doi.org/10.1038/nmeth.3869)), where the user [configuration file](https://github.com/McAllister-NOAA/REVAMP#pipeline-configuration-file--p) defines settings for read filtering based on read length and quality (`filterAndTrim`). After user assessment of success, quality-filtered reads are passed through DADA2, to learn error rates (`learnErrors`), dereplicate sequences (`derepFastq`), merge paired reads (`mergePairs`), check for chimeras (`removeBimeraDenovo`), and assign unique amplicon sequence variants (ASVs), writing out a fasta file and counts feature table. 

Next, ASVs are blasted against the NCBI `nt` database (BLASTn; `subject_besthit`; `max_target_seqs 4000`) ([Camacho et al., 2009](https://doi.org/10.1186/1471-2105-10-421)), exporting percent identity, length of hit, subject taxonomy IDs (taxIDs), and subject accession in tab-delimited format. Alternatively, the pipeline can accept ASV BLASTn results from a user’s custom database, as long as output is in an identical format and taxIDs match that of NCBI’s `taxonomy` database. After the tab-delimited BLAST output file is created, or supplied by the user, the file is then reformatted to simplify the results to only include taxIDs from all best percent identity matches longer than a user-supplied length cutoff. The final step in the pipeline before table and figure generation is taxonomy assignment. Taxonomy assignment with REVAMP is intentionally conservative given: 1) the reference databases are incomplete for some taxa and heavily sampled for others, and 2) markers vary in their ability to resolve different taxa to species (e.g. [Gold et al., 2021](https://doi.org/10.1111/1755-0998.13450)). Before assessing taxonomies, TaxonKit (v 0.5.0) ([Shen and Ren, 2021](https://doi.org/10.1016/j.jgg.2021.03.006)) is run on all discovered taxIDs in the reformatted BLASTn file, with the results reformatted to include kingdom, phylum, class, order, family, genus, and species (`K/P/C/O/F/G/S`) assignments only. In cases where an intermediate taxon is missing in this string, one option is to automatically fill from a lower taxonomy assignment: for example, a `K/P/-/O/F/G/S` string would be replaced with `K/P/O__c/O/F/G/S` to indicate that the class of interest contains the order below it. In REVAMP, taxonomy is then assigned by merging the taxonomic hierarchy of all best BLASTn taxID hits to the lowest common ancestor (i.e. deepest identical taxonomic assignment). This means, for example, that if the best BLASTn hits match several species from the same genus, the ASV will only be identified to the genus level. The only other factor influencing the assignment is the confidence of the taxonomic depth as determined by percent identity cutoffs for each taxonomic level. While different clades vary widely in this respect ([Gebhardt and Knebelsberger, 2015](https://doi.org/10.1007/s10152-015-0434-7); [Zhang and Bu, 2022](https://doi.org/10.3390/insects13050425)), users can choose a set of taxonomic cutoffs per marker gene (see [below](https://github.com/McAllister-NOAA/REVAMP#pipeline-configuration-file--p) for suggestions). For example, if the percent identity falls below a genus assignment cutoff and above the family level cutoff, the final taxonomic assignment can only be to the family level. If SILVAngs ([Glöckner et al., 2017](https://doi.org/10.1016/j.jbiotec.2017.06.1198)) is run on the ASVs file, taxonomic assignments can be made by importing these results, where SILVA assignments are accepted for Bacteria, Archaea, and Eukarya (NCBI taxonomy used for Chloroplast and Mitochondrial hits).

Data exploration and visualization are an integral part of the REVAMP pipeline, as these play an important role in pattern observation and hypothesis generation for follow-up testing. A complete list of products can be found [here](https://github.com/McAllister-NOAA/REVAMP#revamp-results). In general and where appropriate, all figures are generated independently (in separate folders) for comparison of raw read count vs. relative abundance and for ASV-based biological units vs. biological units merged on identical taxonomic assignments. The pipeline produces the following types of figures: 1) [KRONA plots](https://github.com/McAllister-NOAA/REVAMP#00_krona_plots) - an interactive hierarchical data browser ([Ondov et al., 2011](https://doi.org/10.1186/1471-2105-12-385)); 2) [Maps](https://github.com/McAllister-NOAA/REVAMP#01_maps) - simple, bathymetric, and pie-chart inclusive maps are generated using R packages mapdata, marmap, and scatterpie; 3) [Barcharts](https://github.com/McAllister-NOAA/REVAMP#02_barcharts); 4) [Heatmaps](https://github.com/McAllister-NOAA/REVAMP#03_heatmaps); 5) [Alpha diversity metrics](https://github.com/McAllister-NOAA/REVAMP#04_alpha_diversity) - tables and figures sorted and colored by metadata columns; 6) [Ordination plots](https://github.com/McAllister-NOAA/REVAMP#05_ordination) - includes both non-metric multidimensional scaling (NMDS) and principal coordinate analysis (PCoA) plots encircled using sites/replicates/grouping information and colored by environmental metadata; also produces a biplot of samples and taxonomy-colored ASVs and a version faceted by phylum; 7) [Networks](https://github.com/McAllister-NOAA/REVAMP#06_network) - Bray-Curtis distance networks at 0.1 unit increments showing the connectivity patterns of samples and ASVs; 8) [Rarefaction curves](https://github.com/McAllister-NOAA/REVAMP#07_rarefaction_curves) - to estimate sample diversity and effectiveness of sequencing depth; (3-8 are produced through the R package phyloseq \[[McMurdie and Holmes, 2013](https://doi.org/10.1371/journal.pone.0061217)\]); 9) [Environmental fit ordination](https://github.com/McAllister-NOAA/REVAMP#08_environmentfit_ordination) - multivariate statistical analysis of continuous environmental and ASV relative abundance data, plotted in ordination space by the R package vegan ([Oksanen et al., 2020](https://cran.r-project.org/web/packages/vegan/)). Many auxiliary R packages make these figures possible, see [GitHub documentation](https://github.com/McAllister-NOAA/REVAMP#install-the-necessary-dependencies-for-revamp) for more information.

Each stage in read processing is meant to be user guided, since default settings may not apply to each data set. However, if a particular configuration works for the user, all user checkpoints can be avoided (`-y` option). Additionally, should changes need to be made during the pipeline (or after it is run), the user needs to simply revert the progress.txt file to the last usable checkpoint and make any necessary changes before rerunning. Due to REVAMP’s modular structure, the entire pipeline does not need to be rerun when updating taxonomic assignments with an improved reference database.

In addition to the primary workflow of REVAMP, many optional paths are available, including: 1) when provided, information on samples from the same sites or replicates is used to produce figures/statistics averaged over those samples. For [replicate samples](https://github.com/McAllister-NOAA/REVAMP#readsvsreplicatedetection-if-replicates-metadata-column-given), presence/absence tables are generated that can be imported into occupancy modelers (e.g. [unmarked](https://doi.org/10.18637/jss.v043.i10)) and violin plots assessing replicate detection and read abundance are generated; 2) information on positive and negative controls is used to produce visualizations; 3) a negative ASV list can be provided to remove spurious/contaminating ASVs; 4) a positive list (from one taxonomic level) can be provided to produce figures focused on one set of [taxonomies of interest](https://github.com/McAllister-NOAA/REVAMP#taxa_of_interest-if-set-in-figure-configuration-file--f) (e.g. all orders of Copepoda); 5) Taxonomy values of “NA” can be ignored in all figure output. These are in addition to the use of SILVAngs and a user’s custom BLASTn database, as previously discussed above. Together, these options give users simple controls over figure generation.

REVAMP is packaged with four [stand-alone scripts](https://github.com/McAllister-NOAA/REVAMP#stand-alone-applications) that add functionality outside of the main pipeline. These include scripts that allow for 1) independent incorporation of SILVAngs output into the pipeline for the production of tables and figures only (excludes initial ASV assessment); 2) merging SILVAngs Bacterial/Archaeal assignments with BLASTn-based Eukaryotic assignments from two runs of the pipeline for a more complete assessment of all three lineages; 3) production of REVAMP output files from a simple imported table with taxonomy assignments (e.g. binomial species name) and feature counts data (e.g. biomass or density measurements) per sample, such as what could be produced from a traditional morphology-based environmental assessment. Each of these taxonomic assignments is traced to NCBI’s `taxonomy` database through TaxonKit (`name2taxid`), with user guidance; 4) comparison of marker genes (including the output of stand-alone script #3) by merging datasets and creating new tables and figures. These include: 4A) tables comparing counts and identity for shared taxonomies between markers; 4B) Venn diagrams to compare taxonomy counts between markers; 4C) merged taxonomy-based Bray-Curtis distance network figures comparing samples; 4D) merged taxonomy-based ordination plots. All of these products allow the user to assess marker similarity and overlap/novelty.

## Installation
### Easy Installation
Docker container in development

### Installation outside of Docker
REVAMP was developed and tested on MacOS 12.6.3. Additional portability is under development with a Docker container.

```
cd (location of choice: e.g. /Users/mcallister/software/)
git clone https://github.com/McAllister-NOAA/REVAMP.git
export PATH=/Users/mcallister/software/REVAMP:$PATH
```

Confirm that `revamp.sh` is in your `$PATH`: `which revamp.sh`.

#### Install NCBI *nt* and *taxonomy* databases
The blastn step of the pipeline can be completed on another computer, such as a high-performance computing resource. In this case, ASVs.fa can be moved to that resource and `blastn` run there. In either case, the NCBI *nt* and *taxonomy* databases need to be available on whatever resource is processing the `blastn` step.

NCBI databases are prepped for use in REVAMP using the `ncbi_db_cleanup.sh` script.

Dependencies for script: [blast+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/), wget, and [taxonomy-tools](https://github.com/pmenzel/taxonomy-tools), all in your `$PATH`.

Run `ncbi_db_cleanup.sh` from inside blastdb directory of choice:
```
mkdir /path/to/blastdb
cd /path/to/blastdb
ncbi_db_cleanup.sh
```

If NCBI *nt* database files continually fail to download, try to run `update_blastdb.pl nt` outside of the pipeline's script.

#### Install the necessary dependencies for REVAMP

NCBI databases (should be already installed; see above):
* *nt* database (download through `update_blastdb.pl` from the [BLAST+ executables](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)) 
* *taxonomy* database (`wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz`)

Stand-alone tools:
* wget ([HomeBrew](https://formulae.brew.sh/formula/wget) or [MacPorts](https://ports.macports.org/port/wget/) on Mac)
* R (Rscript) - v.4.0.3
* [Taxonomy Tools](https://github.com/pmenzel/taxonomy-tools)
* [Cutadapt](https://cutadapt.readthedocs.io/en/stable/installation.html) - v.3.7
* [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) - v.2.13.0+
* [TaxonKit](https://bioinf.shenwei.me/taxonkit/download/) - v.0.5.0
* [Krona](https://github.com/marbl/Krona/wiki) - v.2.7.1

R packages:
* [dada2](https://benjjneb.github.io/dada2/dada-installation.html) – v.1.18.0
* [phyloseq](https://joey711.github.io/phyloseq/install.html) - v.1.34.0
* vegan - v.2.5-7
* spatstat - v.1.64-1
* dplyr - v.1.0.3
* ggplot2 - v.3.3.3
* ggrepel - v.0.8.2
* ggpubr - v.0.4.0
* ggalt - v.0.4.0
* mapdata - v.2.3.0
* [marmap](https://epante.wordpress.com/marmap/) - v.1.0.5
* mapproj - v.1.2.7
* scatterpie - v.0.1.5
* janitor - v.2.1.0
* VennDiagram - v.1.7.1

Perl dependencies:
* [List::Util](https://metacpan.org/pod/List::Util)
* [List::MoreUtils](https://metacpan.org/pod/List::MoreUtils)

## Quick Start Guide

### Arguments
```
Usage: revamp.sh
       -p Config File
       -f Figure config file
       -s Sample metadata file
       -r Read folder
       -o Output directory
       -b User-supplied BLASTn btab result file (optional)
       -e Toggle use of SILVAngs taxonomy assignments by ASV (optional)
       -y Bypass all terminal prompts (optional)
       -k Remove all intermediate files (optional)
          Not recommended for first run
          (best to rerun with -k after successful completion)
```

Description of arguments:
* `-p`: Pipeline configuration file, including controls for front end quality control, ASV generation, and taxonomy assignment.
* `-f`: Figure configuration file, including controls for filtering, NA display, and identifying taxa of interest.
* `-s`: Sample metadata file, linking sample names with location, grouping information, and metadata (e.g. chemical/physical discreet or continuous parameters).
* `-r`: Folder containing forward (\*_R1.fastq.gz) and reverse (\*_R2.fastq.gz) reads matching sample names in `-s`.
* `-o`: Folder used for all pipeline outputs and organization.
* `-b`: **OPTIONAL** Tab-delimited BLASTn file generated by user. Should be run on the ASVs assigned by DADA2 for the current pipeline run.
* `-e`: **OPTIONAL** Toggle import of SILVAngs results. Should be run on the ASVs assigned by DADA2 for the current pipeline run.
* `-y`: **OPTIONAL** Bypassing all terminal prompts is recommended if the dataset has been run before and chosen parameters work well.
* `-k`: **OPTIONAL** Run REVAMP with this parameter after successful completion to free up disk space by removing intermediate files.

### Required input files

#### Pipeline Configuration File (`-p`)
Modify example file to fit your needs. Primers can contain any base designation from the [IUPAC nucleotide code](https://www.bioinformatics.org/sms/iupac.html) in addition to I (inosine), which is considered as an "N" position. The recommendation for the `blastLengthCutoff` is 90% of the total marker target length. Suggested default `taxonomyCutoffs` are (97,95,90,80,70,60) for rRNA genes and (95,92,87,77,67,60) for protein encoding genes. The `blastMode` parameter has three options, with `allEnvOUT` and `mostEnvOUT` using negative taxID exclusion lists in the BLASTn step. `mostEnvOUT` removes subjects with "uncultured", "environmental samples", "metagenome", and "unidentified" in their "scientific name" field. `allEnvOUT` additionally removes subjects with "unclassified" in the "scientific name" field. Optional parameter `failedMerge_useDirection` should only be used if a prior run of the pipeline has shown that the forward/reverse reads will not merge properly. Optional parameter `removeASVsFILE` can be used to remove specified ASVs (one per line) from analysis (useful for known contaminants). Controls on DADA2 parameters are provided.

```
##This is the REVAMP config file

####################################################
##Frequently modified parameters
####################################################
primerF=GYGGTGCATGGCCGTTSKTRGTT
primerR=GTGTGYACAAAGGBCAGGGAC
blastLengthCutoff=370 #bp length under which BLAST hits are not considered
systemmemoryMB=32000 #MB of total system memory; uses only 70% of max
locationNTdatabase=/path/to/NCBI/blastdb
taxonomyCutoffs=97,95,90,80,70,60 #Percent ID cutoffs for ID to S,G,F,O,C,P
failedMerge_useDirection=FORWARD #Use only FORWARD (R1) or REVERSE (R2) reads (OPTIONAL)
removeASVsFILE=/path/to/file/remove_contaminants.txt (OPTIONAL)

###DADA2 Filtering options:
dada_minlength=120
dada_phix=TRUE
dada_trunQ=2
dada_maxEE1=5
dada_maxEE2=5
dada_trimRight=0 #Recommended to look at fastq quality and trim ends of sequence data as needed.
dada_trimLeft=0
###

####################################################
##Infrequently changed default parameters
####################################################
blastMode=mostEnvOUT #options: allIN, allEnvOUT, mostEnvOUT
```

#### Figure Configuration File (`-f`)
Modify example file to fit your needs. `filterPercent` controls the cutoff relative abundance that is used to place taxa into "zzOther" to simplify bar chart figures. `pieScale` is an attempt to provide some control over the size of pie charts in the "mapPies"-labelled maps. `removeNA` removes any taxonomic assignment of "NA" from bar chart figures. If `filterLowQualSamples` is "TRUE", then `filterPercentLowQualSamples` sets the cutoff for removing a sample from analysis, based on the average percentage of total reads of all other samples. Defining taxa of interest (`providedTaxaOfInterest` "TRUE") can provide many additional figures focusing on only the listed taxa. The `taxaOfInterestFile` must contain a list of the taxa of interest (one per line), all from the same taxonomic level (in this case designated as "Order" for `taxaOfInterestLevel`).

```
##This is the REVAMP figure config file
###Figure options:
filterPercent=5
pieScale=5.7
removeNA=FALSE

###
filterLowQualSamples=TRUE
filterPercentLowQualSamples=10

#TAXA of Interest: Provides specific figures.
#List one name per line. All from same taxonomic level.
providedTaxaOfInterest=TRUE
taxaOfInterestFile=/path/to/file/choice_taxa.txt
taxaOfInterestLevel=Order #Options: Kingdom, Phylum, Class, Order, Family, Genus, or Species
```

#### Sample Metadata File (`-s`)
This important file links a sample to its metadata. At a minimum, users must provide `Sample` designations, with `lat` and `long` recommended. `Sample` should be the first column; however, the order of other columns does not matter. Sample order (i.e. row order) is used for organizing figures. In addition to the controlled vocabulary (below), a user can provide an unlimited number of groupings (applied to visualizations, including NMDS/PCoA encircling) using the column headers from `group1` to `group#`. Any columns that do not match the controlled vocabulary or `group#` are assigned to metadata and can be discrete or continuous variables (e.g. chemical/physical measurements, other observations). Any missing data should be filled in as "NA". Note that while phyloseq can handle missing data, vegan does not calculate environmental fitting with missing data. We recommend that the pipeline is run first with all samples/data columns, followed up with a run to test environmental fitting with the samples with missing metadata and/or metadata columns with missing data removed.

Controlled vocabulary:
* `Sample` = Sample name identical to the "\*" in `*_R1.fastq.gz` and `*_R2.fastq.gz` as found in the reads directory (`-r`). All sample names have "MP_" appended to the front and have non-alphanumeric characters converted to underscores to avoid programmatic errors (e.g. R does not accept row/column headers that start with numbers).
* `lat` = latitude decimal degrees
* `long` = longitude decimal degrees
* `replicates` = String designating replicate samples. Replicate samples are treated both independently and as an averaged aggregate for some figures. Additional replicate-specific figures and tables are created, in particular for occupancy modelling outside the pipeline.
* `sites` = String designating samples taken from the same sites. Similar treatment to the `replicates` column in relation to averaged aggregate visualization.
* `controls` = Can be `negative`, `positive`, or `NA`

User-supplied:
```
Sample	lat	long	replicates    sites  controls      group1 group2 pH_avg temp_degC     water_col_location
18cp45	47.594567	-124.497067	CL42   CL42	NA     AK_trip1      1      7.4    11.4   surface
18cp46	47.594567	-124.497067	CL42   CL42	NA     AK_trip1      1      7.4    11.4   surface
18cp6	47.594567	-124.497067	CL42B  CL42	NA     AK_trip1      2      7.15   3.2    deep
E25_2B_DY20	59.90486667	-171.6984167	M5     M5     NA     AK_trip2      3      7.2    15.1   surface
5-18S_S5_L001	33.97368541	-118.9126162	B5     B5     NA     CA     3      7.22   21.2   deep
positive_control     NA     NA     NA     NA     positive      NA     NA     NA     NA     NA
negative_control     NA     NA     NA     NA     negative      NA     NA     NA     NA     NA
```

Will be modified to:
```
Sample	lat	long	replicates    sites  controls      group1 group2 pH_avg temp_degC     water_col_location
MP_18cp45	47.594567	-124.497067	CL42   CL42	NA     AK_trip1      1      7.4    11.4   surface
MP_18cp46	47.594567	-124.497067	CL42   CL42	NA     AK_trip1      1      7.4    11.4   surface
MP_18cp6	47.594567	-124.497067	CL42B  CL42	NA     AK_trip1      2      7.15   3.2    deep
MP_E25_2B_DY20	59.90486667	-171.6984167	M5     M5     NA     AK_trip2      3      7.2    15.1   surface
MP_5_18S_S5_L001	33.97368541	-118.9126162	B5     B5     NA     CA     3      7.22   21.2   deep
MP_positive_control     NA     NA     NA     NA     positive      NA     NA     NA     NA     NA
MP_negative_control     NA     NA     NA     NA     negative      NA     NA     NA     NA     NA
```

#### Reads directory (`-r`)
Contains forward and reverse reads in the format `*_R1.fastq.gz` and `*_R2.fastq.gz`, respectively, where "\*" matches sample names in the sample metadata (`-s`) file. Fastq files must be gzipped.

### Optional input files

#### User-generated tab-delimited BLASTn file (`-b`)
After ASV-generation, a user can stop the pipeline to run their own BLASTn (use `-outfmt '6 qseqid pident length staxids sacc'`) file on a custom database. To restart, delete the entries in `progress.txt`, leaving only `cutadaptFinished=TRUE` and `dada2_Finished=TRUE` and restart the pipeline.

#### Import of SILVAngs results (`-e`)
This flag will trigger two prompts for file locations:
* `"Enter the location of the SILVAngs ssu or lsu results directory (i.e. ~/Downloads/results/ssu)"`. This folder will contain the SILVAngs export file `~/Downloads/results/ssu/exports/*---otus.csv` and the SILVAngs cluster file `~/Downloads/results/ssu/stats/sequence_cluster_map/data/*.fa.clstr`.
* `"Enter the location of the reference taxonomy map for current SILVA database: i.e. tax_slv_ssu_138.1.txt"`. This file can be downloaded from the Arb SILVA [archive](https://www.arb-silva.de/no_cache/download/archive/current/Exports/taxonomy/).

## REVAMP Results

### REVAMP Directory Structure
```
OUTDIR
├── ASV2Taxonomy
├── blast_results
├── cutadapt
├── dada2
├── processed_tables
│   └── replicate_based_detection (if replicates given)
├── run_logs
└── Figures
    ├── 00_KRONA_plots
    ├── 01_Maps
    ├── 02_Barcharts
    │   └── read_count & relative_abundance
    ├── 03_Heatmaps
    │   └── ASV_based & Taxonomy_merge_based
    ├── 04_Alpha_Diversity
    │   └── ASV_based & Taxonomy_merge_based
    ├── 05_Ordination
    │   ├── ASV_based: read_count & relative_abundance
    │   └── Taxonomy_merge_based
    │       ├── filterInclude_TOSPECIES_only
    │       │   └── read_count & relative_abundance
    │       └── read_count & relative_abundance
    ├── 06_Network
    │   ├── ASV_based: read_count & relative_abundance
    │   └── Taxonomy_merge_based: read_count & relative_abundance
    ├── 07_Rarefaction_Curves
    ├── 08_EnvironmentFit_Ordination
    │   └── ASV_based & Taxonomy_merge_based
    ├── ReadsVSReplicateDetection (if replicates given)
    └── Taxa_of_interest (if given)
        ├── 02_Barcharts
        │   └── read_count & relative_abundance
        ├── 03_Heatmaps
        │   └── ASV_based & Taxonomy_merge_based
        └── 06_Network
            ├── ASV_based: read_count & relative_abundance
            └── Taxonomy_merge_based: read_count & relative_abundance
```

All information printed to screen during a run is stored to `run.log`. If more than one run is done on the same folder, the directory `run_logs` is created and the previous `run.log` is renamed to append the date/time and moved to that directory.

One feature of the REVAMP pipeline is the ability to continue from checkpoints, which is why data files produced by different checkpoints are maintained in their individual folders. Completion of each checkpoint is indicated in the `progress.txt` file, which sets environmental variables for REVAMP.

Checkpoints (add or delete checkpoints to bypass on subsequent pipeline runs):
* `cutadaptFinished=TRUE`: Raw reads are quality controlled using `cutadapt` to remove adaptors and primers (`cutadapt` directory) 
* `dada2_Finished=TRUE`: Trimmed reads are analyzed by `dada2`, where they are quality trimmed, filtered, dereplicated, and merged before ASV assignment (`dada2` directory)
* `blastFinished=TRUE`: ASVs are then "blasted" using BLASTn against NCBI's *nt* (`blast_results` directory) 
* `blastformattingFinished=TRUE`: The btab out file is reformatted to assign best blast hits to each ASV (`blast_results` directory)
* `taxonomyscriptFinished=TRUE`: ASVs are assigned to taxonomy (`ASV2Taxonomy` directory)
* `figuresFinished=TRUE`: Various table products (`processed_tables` directory) and figures (`Figures` directory) are generated.
Note: If you wish to redo a checkpoint that has already passed, simply delete sequentially from the bottom of the `progress.txt` file until you have deleted the step you want to redo. **It is not recommended to delete a step in the middle while allowing the other steps to remain. It will break.**

The `Figure` directory includes four primary choices for preference of data analysis (where applicable). The User can choose to look at figures where ASVs are considered a single ecological unit (`ASV_based`) or where ASVs have been collapsed so that each represents a single unique taxonomic hierarchy (`Taxonomy_merge_based`). Further, figures can be generated based on the raw count data (`read_count`) or based on the normalized relative percent abundance (`relative_abundance`).

### REVAMP Files of Interest

All configuration files are copied to the REVAMP results folder. `config_file.txt` controls all of the main options at the front end of REVAMP (is a copy of the file provided by `-p`. `figure_config_file.txt` controls all the options for the back-end figure generation of REVAMP (is a copy of the file provided by `-f`. The provided sample metadata file is copied (`sample_metadata.txt`) and converted/cleaned for use in R (`sample_metadata_forR.txt`).

Each R script outputs stdout and stderror to a log file within the R scripts primary directory. In addition, if you wish to customize or debug any R-based figure/file, you can open any of the R scripts in the `assets` folder. The second field of the tab-delimited file `Rscript_arguments.log` gives each of the positional arguments fed to each script during the pipeline. Simply uncomment the args block at the front end of the script and input each argument.

Besides the outputs in the `processed_tables` and `Figures` directory, the User may also find these files useful:

* `OUTDIR/sample_order.txt` - This file is created from the sample metadata file (`-s`), though it can be modified if the user desires a different order for figures
* `OUTDIR/dada2/ASVs.fa` - This is the fasta file with ASV nucleotide sequences used in all subsequent steps
* `OUTDIR/dada2/ASVs_counts.tsv` - This biom table is the original table showing read counts per ASV per sample
* `OUTDIR/blast_results/ASV_blastn_nt.btab` - BLASTn tab-delimited output file
* `OUTDIR/blast_results/ASV_blastn_nt_formatted.txt` - Shows the list of identical best hits for each ASV, with the taxIDs listed
* `OUTDIR/ASV2Taxonomy/basic_ASV_taxonomy_stats.txt` - Counts for the number of ASVs by sample and the numbers of ASVs assigned to different taxonomic levels
* `OUTDIR/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy.tsv` - This biom table shows the read counts per ASV per sample, where unique taxonomy is collapsed (summed) to be represented by only a single ASV
* The ASV counts files are also represented in `OUTDIR/ASV2Taxonomy` with "Unknowns" removed
* `OUTDIR/ASV2Taxonomy/outname_asvTaxonomyTable.txt` - This is the primary record for taxonomy assignments for each ASV (also available with "Unknowns" removed
* `OUTDIR/ASV2Taxonomy/reformatted_taxonkit_out.txt` - This file shows all taxonomy records for taxIDs in the reformatted BLASTn file. Can be manually modified in the pipeline and/or gaps automatically filled as described.
* `OUTDIR/ASV2Taxonomy/taxonomy2PercentAbundance_humanReadable.txt` - Relative abundance (%) for each unique taxonomic string by sample. In this file, numbers are rounded to two decimal places and any value less than 0.01% rounded is marked with a "D" to indicate detection.
* `OUTDIR/ASV2Taxonomy/taxonomy2PercentAbundance_humanReadable_NoRounding.txt` - Same as the file above, but values are not rounded (and there are no "D" designations)
* `OUTDIR/ASV2Taxonomy/outname_singleBlastHits_with_MULTItaxid.txt` - File showing the ASVs where at least one of the best blast hit sequences is marked as belonging to two or more taxonomic IDs in NCBI. REVAMP will consider the last common ancestor of such a hit if it is the only available (marked as "TRUE" in the "USED" column). However, if other single taxID hits are available, REVAMP with use those assignments instead (marked as "FALSE" in the "USED" column).
* `OUTDIR/ASV2Taxonomy/OUTDIR_heatmap_multiASV.txt` - Shows a heatmap of read counts for taxonomic assignments with more than one ASV. This kind of result could show where multi-ASVs with the same taxonomy might be the result of erroneous sequence/chimeric error (i.e. one ASV with lots of reads with a few others with the same taxonomy but very few sporadic read hits) as well as which multi-ASVs might represent different ecotypes of the same taxa (i.e. multiple ASVs with believable read distribution among samples).
* `OUTDIR/ASV2Taxonomy/OUTDIR_unknown_asvids.txt` - File can help to track the reasons for "Unknown" assignments

#### REVAMP `processed_tables` Folder

Any counts file modified after initial assignments in DADA2 is stored in the `processed_tables` folder. The files created depend on whether sample metadata columns for `controls`, `replicates`, and/or `sites` are provided.

At a minimum:
* `ASVs_counts_NOUNKNOWNS_percentabund.tsv` = Read counts converted into relative percent abundance with ASVs with unknown taxonomic assignments removed from total.
* `ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.tsv` = Same as above, except that reads across ASVs with identical taxonomic assignments are summed prior to the calculation. One ASV name is chosen as the representative for the group, and the rest are removed from the file.
* `ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt` = ASV taxonomy table built to simplify figures based on `filterPercent`, where the taxonomy at any hierarchical level is assigned to "zzOther" if it falls below the `filterPercent` cut off. 

If a `controls` metadata column is given, those samples are removed from downstream analysis:
* `ASVs_counts_controlsRemoved.tsv` = ASV counts with the control samples removed
* `ASVs_counts_NOUNKNOWNS_controlsRemoved.tsv` = ASV counts with control samples and unknown ASVs removed
* `ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_controlsRemoved.tsv` = ASV counts collapsed on identical taxonomic assignments, with control samples removed

If either `replicates` or `sites` metadata columns are given, then an additional set of files is created, whereby samples with identical `replicates` or `sites` labels are grouped together by averaging relative abundance. This also concatenates the sample metadata file, combining and dereplicating entries for each metadata column.

If either `replicates` or `sites` metadata columns are given:
* `ASVs_counts_NOUNKNOWNS_percentabund_groupedByReplicates.tsv`
* `ASVs_counts_NOUNKNOWNS_percentabund_groupedBySites.tsv`
* `ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund_groupedByReplicates.tsv`
* `ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund_groupedBySites.tsv`
* `sample_metadata_NOUNKNOWNS_percentabund_groupedByReplicates.tsv`
* `sample_metadata_NOUNKNOWNS_percentabund_groupedBySites.tsv`

In addition, the following files are created in the `replicate_based_detection` folder if a `replicates` sample metadata column is given:
* `compRelAbund_replicateDetection_ASVbased_NoUnknowns_filtsamples.txt` = Stats file on each ASV counts by replicates without ASVs with an unknown taxonomic assignment.
* `compRelAbund_replicateDetection_ASVbased_withUnknowns_filtsamples.txt` = Stats file on each ASV counts by replicates with ASVs with an unknown taxonomic assignment.
* `compRelAbund_replicateDetection_TAXAbased_NoUnknowns_filtsamples.txt` = Stats file on ASV counts by replicates for ASVs collapsed on identical taxonomic assignments, no unknowns.
* `presenceabsence_unmarked_ASVbased_NoUnknowns_filtsamples.txt` = Input file for Unmarked with presence/absence indicated for each ASV by replicate. No unknowns.
* `presenceabsence_unmarked_ASVbased_withUnknowns_filtsamples.txt` = Input file for Unmarked with presence/absence indicated for each ASV by replicate. With unknowns.
* `presenceabsence_unmarked_TAXAbased_NoUnknowns_filtsamples.txt` = Input file for Unmarked with presence/absence indicated for each ASV by replicate. No unknowns. Collapsed on taxonomy.

#### REVAMP `Figures` Folder

```
└── Figures
    ├── 00_KRONA_plots
    ├── 01_Maps
    ├── 02_Barcharts
    │   └── read_count & relative_abundance
    ├── 03_Heatmaps
    │   └── ASV_based & Taxonomy_merge_based
    ├── 04_Alpha_Diversity
    │   └── ASV_based & Taxonomy_merge_based
    ├── 05_Ordination
    │   ├── ASV_based: read_count & relative_abundance
    │   └── Taxonomy_merge_based
    │       ├── filterInclude_TOSPECIES_only
    │       │   └── read_count & relative_abundance
    │       └── read_count & relative_abundance
    ├── 06_Network
    │   ├── ASV_based: read_count & relative_abundance
    │   └── Taxonomy_merge_based: read_count & relative_abundance
    ├── 07_Rarefaction_Curves
    ├── 08_EnvironmentFit_Ordination
    │   └── ASV_based & Taxonomy_merge_based
    ├── ReadsVSReplicateDetection (if replicates given)
    └── Taxa_of_interest (if given)
        ├── 02_Barcharts
        │   └── read_count & relative_abundance
        ├── 03_Heatmaps
        │   └── ASV_based & Taxonomy_merge_based
        └── 06_Network
            ├── ASV_based: read_count & relative_abundance
            └── Taxonomy_merge_based: read_count & relative_abundance
```

##### 00_KRONA_plots
<img width="1741" alt="Krona Plot" src="https://github.com/McAllister-NOAA/REVAMP/assets/60410177/2b17c26c-337d-4ccf-8cc7-9cbd0b968170">

Interactive hierarchical data browser allows users to explore the biodiversity of each sample (`out_master_krona.html`) or summed reads from all samples in the run (`out_samplesSummedKRONA.html`) at any taxonomic level they desire. In the upper right corner, a summary of the read count for each selected taxa is displayed. Results are searchable.

##### 01_Maps
![maos](https://github.com/McAllister-NOAA/REVAMP/assets/60410177/fa2477b1-d84e-4d34-94d5-15659924d038)

Three types of maps are created from the latitude and longitude data given in the sample metadata file, including a basic map, a bathymetric map, both with data points (scattered for visibility), and a basic map with the data points replaced with pie charts showing the unique terminal taxa found in each sample.

##### 02_Barcharts
Bar charts are useful for exploring both read depth (sequencing effort) and patterns in relative abundance. Note that relative sequence abundance does not necessarily equate to relative abundance of an organism in the sample, as several factors can influence read depth.

Many bar charts are produced by the pipeline, separated into read-count-centric and relative-abundance-centric figures. Figures are created at different taxonomic levels (Phylum to Species), as well as for all unique terminal taxa (no matter the depth). In addition, the `filterPercent` from the figure configuration file (`-f`) is used to simplify figures by placing taxa less than the designated percentage of total community into the "zzOther" category. This can make viewing of the figure much easier. If there are still too many colors to ascertain the differences between taxa, one method for exploring the bar charts is to load them (and their legend) into Adobe Illustrator and to use the "Select"/"Same"/"Fill Color" feature to select all the same taxa.

|![barcharts-01](https://github.com/McAllister-NOAA/REVAMP/assets/60410177/a8a09488-b6aa-4c35-a3ed-c985622ebfa1)|
|:--:|
|Bar charts comparing read-count- (`barplot_readcount_allsamples_alltaxa_Class.pdf`) and relative-abundance-centric (`barplot_relabund_filtsamples_alltaxa_Class.pdf`) figures|

|![barcharts-03](https://github.com/McAllister-NOAA/REVAMP/assets/60410177/7569b6b7-7da3-421b-9067-4ea66924a612)|
|:--:|
|Bar chart showing site-averaged relative abundance (`barplot_relabund_siteGroupedSamples_alltaxa_Class.pdf`)|

|![barcharts-02](https://github.com/McAllister-NOAA/REVAMP/assets/60410177/7bcc1860-c3e6-446e-a982-e2628f98e928)|
|:--:|
|Bar chart showing unique terminal taxa (`barplot_relabund_allsamples_filtLowAbundTaxa_to_zzOther_uniqueTerminalTaxa_noUnknowns.pdf`). This figure filters taxa less than the designated percent relative abundance to "zzOther", though the legend is shown to highlight the complexity of such a figure if this filtering is not done (`barplot_relabund_allsamples_uniqueTerminalTaxa_noUnknowns_legend.pdf`).|

##### 03_Heatmaps
![heatmap](https://github.com/McAllister-NOAA/REVAMP/assets/60410177/6222b512-cf32-4683-87be-862a3278411c)

Heatmaps can be useful for looking at patterns in ASV presence/absence and relative abundance between samples (`heatmap_ASVs_relabund_filtsamples_alltaxa_clustSamples.pdf`). Figures are created that order samples by the sample metadata file (`-s`), or by clustering samples by Jaccard similarity using the NMDS ordination method. In both cases, ASV's are ordered using the NMDS ordination method based on Jaccard distance. Heatmaps are created for ASV-based and merged-Taxonomy-based datasets. In addition (as in the figure above, right side), the site- or replicate-grouped datasets are used for alternative heatmaps (`heatmap_ASVs_relabund_siteGroupedSamples_alltaxa_clustSites.pdf`).

##### 04_Alpha_Diversity
In addition to creating tables for alpha diversity metrics, REVAMP creates several figures to visualize alpha diversity, including in relation to all groupings and other metadata columns. 

|![alphaDiversity_ASVs_normalized_filtsamples_alltaxa_allMeasures_group3Colored](https://github.com/McAllister-NOAA/REVAMP/assets/60410177/f5032842-d114-4e4f-ac2c-8cf60ab8a768)|
|:--:|
|Alpha diversity on a per sample basis, colored by region (`alphaDiversity_ASVs_normalized_filtsamples_alltaxa_allMeasures_group3Colored.pdf`)|

|![alphaDiversity_ASVs_normalized_siteClusteredSamples_alltaxa_allMeasures_group2Colored](https://github.com/McAllister-NOAA/REVAMP/assets/60410177/9fd9d986-8373-44c6-9419-0aa2c2d390e1)|
|:--:|
|Alpha diversity with samples from the same sites displayed on the same x-axis column while still calculated separately (i.e. site-clustered). Colored by depth. (`alphaDiversity_ASVs_normalized_siteClusteredSamples_alltaxa_allMeasures_group2Colored.pdf`)|

|![alphaDiversity_ASVs_normalized_siteGroupedSamples_alltaxa_allMeasures_group2Colored](https://github.com/McAllister-NOAA/REVAMP/assets/60410177/70e4a0f3-8950-43fe-b407-2eb005aefc5c)|
|:--:|
|Alpha diversity with samples from the same sites averaged and calculated as a single aggregated sample (`alphaDiversity_ASVs_normalized_siteGroupedSamples_alltaxa_allMeasures_group2Colored.pdf`). Since most sites included both surface and deep samples, the concatenated sample metadata is not as useful for this example, but this can highlight patterns of missing data (i.e. red and green dots in the figure are missing representatives from the other depth regime).|

In addition to alpha diversity calculations based on the whole dataset, REVAMP also calculates alpha diversity with subsets of the data that include only those ASVs that have a complete taxonomy to the desired depth (made from Phylum to Species), as well as figures based on the alpha diversity metrics calculated for those ASVs complete to the Species level. 

##### 05_Ordination
Ordination is a very powerful tool for using multivariate statistics to calculate and visualize relationships between samples based on sample composition. Here, both non-metric multidimensional scaling (NMDS; preference) and principal coordinate analysis (PCoA) are used to make ordination plots (as well as tables with coordinates and convex hull analyses). Ordination plots are made in four combinations: 1) ASV-based/read count, 2) ASV-based/relative abundance, 3) Taxonomy-merge-based/read count, and 4) Taxonomy-merge-based/relative abundance.

|![ordination-01](https://github.com/McAllister-NOAA/REVAMP/assets/60410177/7ba8e4cb-6b11-481b-a5db-4e0e9a685c73)|
|:--:|
|Basic NMDS plot (ASV-based/read count) with samples labelled (can be removed in Illustrator) (`NMDS_ASVbased_samples_relabund_filtsamples_alltaxa.pdf`)|

|![ordination-02](https://github.com/McAllister-NOAA/REVAMP/assets/60410177/3e4c425f-edea-40e6-874d-72ebd0cc7dce)|
|:--:|
|Examples of coloring (data points) and encircling (outline of identical metadata categories) of NMDS plot. A) Coloring and encircling based on metadata group (`NMDS_ASVbased_samples_relabund_filtsamples_alltaxa_group3Colored_encircleGroup.pdf`). B) Coloring based on metadata group and encircling sites (`NMDS_ASVbased_samples_relabund_filtsamples_alltaxa_group3Colored_encircleSites.pdf`). C) Coloring based on temperature (`NMDS_ASVbased_samples_relabund_filtsamples_alltaxa_chem_temperature_degC_Colored.pdf`).|

|![ordination-03](https://github.com/McAllister-NOAA/REVAMP/assets/60410177/99d5ac0a-5245-417a-9e92-852cf2b28bad)|
|:--:|
|Biplot of samples and ASVs in same ordination space. Default is to color by Phylum. (`NMDS_ASVbased_taxa_relabund_filtsamples_alltaxa_biplot.pdf`)|

|![ordination-04](https://github.com/McAllister-NOAA/REVAMP/assets/60410177/24249da1-2f11-49b9-9611-1755fada24bd)|
|:--:|
|NMDS plot faceted by Phylum and colored by Genus (`NMDS_ASVbased_taxa_relabund_filtsamples_alltaxa_facetPhylum_colorGenus.pdf`)|

##### 06_Network
Network figures are useful for organizing information about sample relatedness. In the figures produced by REVAMP, network figures are created for samples and ASVs at 0.1 Bray-Curtis distance increments. Only samples with a smaller distance than the distance metric cutoff have a line drawn between them. In this way, users can visualize the relationships of samples and ASVs. In the figures below, the figures at 0.2 and 0.9 distance are not very useful, while those in the middle show clusters of related samples/ASVs.

|![network-01](https://github.com/McAllister-NOAA/REVAMP/assets/60410177/6a60bb4d-45e2-4384-8e90-40b3c423de52)|
|:--:|
|Network figure plotting samples and colored by sample site (if given) (`Network_ASVbased_samples_relabund_filtsamples_gte1perctaxa_dist0.2_SampleLabeled_colorSites.pdf`, etc.)|

|![network-02](https://github.com/McAllister-NOAA/REVAMP/assets/60410177/878776f6-897d-45b2-b9c1-3b1b740b9c1d)|
|:--:|
|Network figure plotting ASVs (`Network_ASVbased_taxa_relabund_filtsamples_gte1perctaxa_dist0.2_ASVLabeled.pdf`, etc.)|

##### 07_Rarefaction_Curves
Rarefaction curves can highlight diversity difference between sample as well as samples with inadequate read depth (below `RarefactionCurve_ASVbased_readcount_allsamples_alltaxa.pdf`).

![RarefactionCurve_ASVbased_readcount_allsamples_alltaxa](https://github.com/McAllister-NOAA/REVAMP/assets/60410177/851833c8-c179-4ceb-b56c-a26e65842bc3)

##### 08_EnvironmentFit_Ordination
The R package `vegan` is used to run the environmental fitting, which means that the NMDS is recalculated by `vegan`. New figures and tables are created in this folder to show the new coordinates so as not to confuse assessment.

|![environfit-01](https://github.com/McAllister-NOAA/REVAMP/assets/60410177/d303123c-29b4-4504-b1f9-ed1a6988c1de)|
|:--:|
|`vegan`-based NMDS showing labelled samples (left; `NMDS_vegan_ASVbased_samples_relabund_filtsamples_alltaxa.pdf`) and encircled samples based on site name (right; `NMDS_vegan_ASVbased_samples_relabund_filtsamples_alltaxa_encircleSites.pdf`).|

Each variable's goodness of fit (r^2) and statistical probability (p-value) is stored in a table. ASV-based or taxonomy-merge-based relative abundance is also tested as a variable, so that organisms significantly structuring the ordination space can be visualized. Figures are generated displaying "chemistry" (non-controlled vocabulary columns) and ASV variables ≤0.05 p-value, "chemistry" and ASV variables ≤0.01 p-value, and all given metadata columns.

|![environfit-02](https://github.com/McAllister-NOAA/REVAMP/assets/60410177/8c13d830-086e-4852-baff-36d15871685b)|
|:--:|
|`vegan`-based NMDS overlayed with A) variables significant at ≤0.05 p-value (`NMDS_vegan_ASVbased_samples_relabund_filtsamples_alltaxa_environmentFitVectors_p0.05.pdf`), B) variables significant at ≤0.01 p-value (`NMDS_vegan_ASVbased_samples_relabund_filtsamples_alltaxa_environmentFitVectors_p0.001.pdf`), and C) all "chemistry" variables (`NMDS_vegan_ASVbased_samples_relabund_filtsamples_alltaxa_environmentFitVectors_chemOnly.pdf`)|

##### ReadsVSReplicateDetection (if `replicates` metadata column given)
In addition to the tables created in the `processed_tables` folder, REVAMP created a violin plot examining relative abundance of reads compared to the number of replicates ASVs are found in. In a perfect world, abundant ASVs should be found in every sample, and less abundant ASVs should not. These figures highlight the stochasticity of metabarcoding data. Figures are created with both ASV-based and Taxonomy-merge-based datasets, using A) averaged read counts over all replicates, B) averaged read counts over the number of replicates where the ASV/taxa was found, and C) summed reads over all replicates.

![ViolinBoxPlot_ReadsVSReplicateDetection_ASVbased_withoutUnknowns_avgReadsAllReplicates](https://github.com/McAllister-NOAA/REVAMP/assets/60410177/8436c30d-06a8-436c-a8ae-9861412a835c)
`ViolinBoxPlot_ReadsVSReplicateDetection_ASVbased_withoutUnknowns_avgReadsAllReplicates.pdf`

##### Taxa_of_interest (if set in Figure configuration file `-f`)
REVAMP can create a subset of figures (where appropriate), if the user provides a list of taxonomies of interest (must be same taxonomic level). The user must set the following parameters in the Figure configuration file (`-f`): `providedTaxaOfInterest=TRUE`, `taxaOfInterestLevel` (Options: Kingdom, Phylum, Class, Order, Family, Genus, or Species), and `taxaOfInterestFile` (one taxa per line). Figures created include only: 02_Barcharts, 03_Heatmaps, and 06_Network.

## Stand-alone applications

### silvangs_convertTable2REVAMP.sh
Allows for the independent incorporation of SILVAngs output into the pipeline for the production of tables and figures only. Unlike the `-e` option in the main `revamp.sh` pipeline, this approach excludes initial quality control through ASV assignment. Ideal for use with single read data (e.g. Nanopore reads), producing the same figures and tables as the main pipeline. 

#### Arguments
```
Usage: silvangs_convertTable2REVAMP.sh"
       -i Input SILVAngs exports/x---[ls]su---otus.csv spreadsheet
       -r Reference taxonomy map for current SILVA database: i.e. tax_slv_ssu_138.1.txt
       -s Sample metadata file
       -o Output directory
       -f Filter percent cutoff for assignment to zzOther
       -n Filter NAs from figures (optional)
       -t Taxa of interest file (one per line) (optional)
       -c Taxonomic category (e.g. Order) used in Taxa of interest file (required if -t called)
       -m Merge NCBI Eukaryote taxa assignments with SILVA Bacteria/Archaea assignments (optional)
       -d CD-HIT clustr mapping file if clustering done before hand (optional)
```

Description of arguments:
* `-i`: SILVAngs export file `~/Downloads/results/ssu/exports/*---otus.csv`
* `-r`: Reference taxonomy map for current SILVA database (i.e. [tax_slv_ssu_138.1.txt](https://www.arb-silva.de/no_cache/download/archive/current/Exports/taxonomy/))
* `-s`: Sample metadata file, as above
* `-o`: Folder used for all pipeline outputs and organization
* `-f`: Percentage cut off for assigning low abundance taxa to zzOther
* `-n`: **OPTIONAL** Toggle for filtering "NA" assignments from phyloseq figures.
* `-t`: **OPTIONAL** Toggle for providing a file of taxa of interest for additional figure generation.
* `-c`: **OPTIONAL;REQUIRED if `-t` is provided** Taxonomic level of the taxa provided in `-t`.
* `-m`: **OPTIONAL** Creates a merged result folder where NCBI Eukaryote taxa assignments are applied from the SILVAngs export file when SILVA assignments are to "Mitochondria" or "Chloroplast".
* `-d`: **OPTIONAL** If read files were clustered with CD-HIT before supplying them to SILVAngs, provide that clustr file here to inflate read counts for clusters.

### mergeBLASTandSILVAngs_ASV2Taxonomy.pl
Allows for merging of two independent runs of the pipeline, one with regular REVAMP with BLASTn-based taxonomic assignments and one using SILVAngs taxonomic assignments (with `-e` in the main pipeline). In this case, SILVAngs Bacterial/Archaeal assignments are prioritized (originally developed and more accurate for these lineages), while Eukaryotic assignments are prioritized from the BLASTn-based approach (yields more specific assignments than SILVAngs). Runs must be on identical ASVs, which can be accomplished by running first the default pipeline, then copying the outdirectory and modifying `progress.txt` to start after the DADA2 checkpoint with `-e` flagged. This script effectively is the same as the `taxonomyscriptFinished=TRUE` checkpoint in `progress.txt`, and the pipeline run in the merged directory can be continued from the main REVAMP pipeline after running it. 

#### Arguments
```
-a = ASV counts table (make sure there is text in the upper left)
-b = ASV taxonomy table from BLAST run
-s = ASV taxonomy table from SILVAngs run
-n = Allin Output basename
-o = List of samples (one per line) in the order you want them exported. Must be exact matches to ASV counts table.
     Does not have to include all samples. (optional)
-m = REVAMP directory
-h = This help message
```

Description of arguments:
* `-a`: ASV counts table (`ASVs_counts.tsv`) from either run (should be the same)
* `-b`: ASV taxonomy table from BLASTn run (`outdir_asvTaxonomyTable.txt`)
* `-s`: ASV taxonomy table from SILVAngs run (`outdir_asvTaxonomyTable.txt`)
* `-n`: Out directory name for the merged run.
* `-o`: **OPTIONAL** Sample order file (normally there is one created in the out directory after REVAMP interprets the sample metadata file.
* `-m`: PATH to the REVAMP program directory.

### morphology_convertTable2REVAMP.sh
Allows for the production of figures and tables from simple taxonomy assignments (e.g. binomial species name or higher) with feature counts data (e.g. biomass or density measurements) per sample, such as what could be used for morphology-based environmental assessment. Each of these taxonomic assignments is traced to NCBI’s taxonomy database through TaxonKit (name2taxid), with user guidance where needed.

#### Arguments
```
Usage: morphology_convertTable2REVAMP.sh
       -i Input morphology spreadsheet
       -s Sample metadata file
       -o Output directory
       -f Filter percent cutoff for assignment to zzOther
       -n Filter NAs from figures (optional)
       -t Taxa of interest file (one per line) (optional)
       -c Taxonomic category (e.g. Order) used in Taxa of interest file (required if -t called)
       -y Automate filling in taxonkit output (recommended; optional)
```

Description of arguments:
* `-i`: Input file in the format "Sample\tTaxonomy\tCount1...CountX" with headers. Samples should match the `-s` sample metadata file. Taxonomy can be binomial species names or higher taxonomic levels. Count data can be any abundance measurement (e.g. density or biomass measure for taxa in the sample), and as many measurements as you wish can be provided.
* `-s`: Sample metadata file, as above
* `-o`: Name of the output directory
* `-f`: Percentage cut off for assigning low abundance taxa to zzOther
* `-n`: **OPTIONAL** Toggle for filtering "NA" assignments from phyloseq figures.
* `-t`: **OPTIONAL** Toggle for providing a file of taxa of interest for additional figure generation.
* `-c`: **OPTIONAL;REQUIRED if `-t` is provided** Taxonomic level of the taxa provided in `-t`.
* `-y`: **OPTIONAL;RECOMMENDED** Automatically fills gaps in the taxonkit output as described in the main pipeline. Where a blank "ORDER" with known "FAMILY" would be assigned "FAMILY__o" to designate that it is the order that contains that family.

### compare_markers_TablesFigures.sh
Allows for the comparison of marker genes and the output from `morphology_convertTable2REVAMP.sh` by merging datasets and creating new tables and figures.
Tables and figures include:
1. Tables comparing counts and identity for shared taxonomies between markers
2. Venn diagrams to compare taxonomy counts between markers
3. Merged taxonomy-based Bray-Curtis distance network figures comparing samples
4. Merged taxonomy-based ordination plots. 

#### Arguments
```
Usage: compare_markers_TablesFigures.sh
       -i Input folder with internal folders: ASV_relabund (req), ASV_taxonomy (req), Sample_metadata (req), and Taxa_relabund_Human (opt)
          Each folder should have their corresponding files renamed to Marker.txt (i.e. COI.txt)
       -o Output directory
       -s Sample equivalents file (Marker tab MarkerSample tab ReferenceSample) (optional)
       -r Reference marker file name (Marker.txt) (required if -s is called)
       -t Taxa of interest file (one per line) (optional)
       -c Taxonomic category (e.g. Order) used in Taxa of interest file (required if -t called)
```

Description of arguments:
* `-i`: Input directory containing the four specifically named folders as stated above. In the `ASV_relabund` folder, copy the file `ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.txt` from the `processed_tables` folder for each marker (relabel as `Marker.txt`). In the `ASV_taxonomy` folder, copy the file `outname_asvTaxonomyTable_NOUNKNOWNS.txt` from the `ASV2Taxonomy` folder for each marker. In the `Sample_metadata` folder, copy the file `sample_metadata_forR.txt` from the main pipeline out directory for each marker. In the `Taxa_relabund_Human` folder, copy the file `taxonomy2PercentAbundance_humanReadable.txt` or `taxonomy2PercentAbundance_humanReadable_NoRounding.txt` from the `ASV2Taxonomy` folder for each marker.
* `-o`: Location of the output directory
* `-s`: **OPTIONAL;REQUIRED if `Taxa_relabund_Human` folder given** Tell the program what samples are equivalent between a marker and a "reference" marker chosen for assessment of the taxonomy2PercentAbundance files. Format is tab-delimited, with three columns: 1) marker name, 2) sample names from the marker, and 3) sample names as they are found in the reference marker.
* `-r`: **OPTIONAL;REQUIRED if `-s` is provided** Name of the reference marker
* `-t`: **OPTIONAL** Toggle for providing a file of taxa of interest for additional figure generation.
* `-c`: **OPTIONAL;REQUIRED if `-t` is provided** Taxonomic level of the taxa provided in `-t`.

## Miscellaneous uses

### Filtering BLAST data by release date
One beneficial application of the REVAMP pipeline is the assessment of marker genes and database quality over time. By filtering NCBI's GenBank sequences from the BLASTn results based on when they were released, REVAMP allows a user to determine the effect of adding reference organisms to the database, and how effective that database is at resolving different marker genes to Species, Genus, Family, Order, and Class. See the docs in this repository: [BLAST_dateFiltering](https://github.com/McAllister-NOAA/BLAST_dateFiltering).

#### Legal Disclaimer
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
