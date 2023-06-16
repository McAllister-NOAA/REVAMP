![REVAMP_logo_v3_color-tagline](https://github.com/McAllister-NOAA/REVAMP/assets/60410177/00f0f999-043d-43c9-8bf8-1b8aa86bebca)

# REVAMP: Rapid Exploration and Visualization through an Automated Metabarcoding Pipeline

Created by Sean McAllister, Chris Paight, Emily Norton, and Matt Galaska

**REVAMP** is designed to streamline the processing of metabarcoding data from raw fastq data files to the generation of information and visuals. The purpose is to standardize methods for rapid assessments of molecular ecology datasets, including those for environmental DNA (eDNA) monitoring.

Developed by the **Ocean Molecular Ecology** ([OME](https://www.pmel.noaa.gov/ocean-molecular-ecology/)) Group at the Pacific Marine Environmental Laboratory ([PMEL/NOAA](https://www.pmel.noaa.gov/)) in cooperation with the Cooperative Institute for Climate, Ocean, & Ecosystem Studies ([CICOES/UW](https://cicoes.uw.edu/)).

## Citation

**If you find REVAMP useful in your research, please cite...**

```
IN REVIEW
```

**REVAMP is a wrapper for severl tools, which should also be cited:**

* Cutadapt – [Martin, 2011](https://doi.org/10.14806/ej.17.1.200)
* DADA2 – [Callahan et al., 2016](https://doi.org/10.1038/nmeth.3869) ([useful tutorial](https://benjjneb.github.io/dada2/index.html))
* BLASTn – [Camacho et al., 2009](https://doi.org/10.1186/1471-2105-10-421)
* TaxonKit – [Shen and Ren, 2021](https://doi.org/10.1016/j.jgg.2021.03.006)
* Krona – [Ondov et al., 2011](https://doi.org/10.1186/1471-2105-12-385)
* phyloseq – [McMurdie and Holmes, 2013](https://doi.org/10.1371/journal.pone.0061217) ([useful tutorial](https://joey711.github.io/phyloseq/))
* vegan – [Oksanen et al., 2020](https://cran.r-project.org/web/packages/vegan/)

Many other auxiliary R packages make the figures possible: see dependencies.

## Installation
### Easy Installation
Docker container in development

### Installation outside of Docker
REVAMP was developed and test on MacOS 12.6.3. Additional portability is under development with Docker container.

```
cd (location of choice: e.g. /Users/mcallister/software/)
git clone https://github.com/McAllister-NOAA/REVAMP.git
export PATH=/Users/mcallister/software/REVAMP:$PATH
```

Confirm that ```revamp.sh``` is in your ```$PATH```: ```which revamp.sh```.

#### Install NCBI *nt* and *taxonomy* databases
The blastn step of the pipeline can be completed on another computer, such as a high-performance computing resource. In this case, ASVs.fa can be moved to that resource and ```blastn``` run there. In either case, the NCBI *nt* and *taxonomy* databases need to be available on whatever resource is processing the ```blastn``` step.

NCBI databases are prepped for use in REVAMP using the ```ncbi_db_cleanup.sh``` script.

Dependencies for script: [blast+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/), wget, and [taxonomy-tools](https://github.com/pmenzel/taxonomy-tools), all in your ```$PATH```.

Run ```ncbi_db_cleanup.sh``` from inside blastdb directory of choice:
```
mkdir /path/to/blastdb
cd /path/to/blastdb
ncbi_db_cleanup.sh
```

If NCBI *nt* database files continually fail to download, try to run ```update_blastdb.pl nt``` outside of the pipeline's script.

#### Install the necessary dependencies for REVAMP

NCBI databases (should be already installed; see above):
* *nt* database (download through ```update_blastdb.pl``` from the [BLAST+ executables](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)) 
* *taxonomy* database (```wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz```)

Stand alone tools:
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
* ```-p```: Pipeline configuration file, including controls for front end quality control, ASV generation, and taxonomy assignment.
* ```-f```: Figure configuration file, including controls for filtering, NA display, and identifying taxa of interest.
* ```-s```: Sample metadata file, linking sample names with location, grouping information, and metadata (e.g. chemical/physical discreet or continuous parameters).
* ```-r```: Folder containing forward (\*_R1.fastq.gz) and reverse (\*_R2.fastq.gz) reads matching sample names in ```-s```.
* ```-o```: Folder used for all pipeline outputs and organization.
* ```-b```: **OPTIONAL** Tab-delimited BLASTn file generated by user. Should be run on the ASVs assigned by DADA2 for the current pipeline run.
* ```-e```: **OPTIONAL** Toggle import of SILVAngs results. Should be run on the ASVs assigned by DADA2 for the current pipeline run.
* ```-y```: **OPTIONAL** Bypassing all terminal prompts is recommended if the dataset has been run before and chosen parameters work well.
* ```-k```: **OPTIONAL** Run REVAMP with this parameter after successful completion to free up disk space by removing intermediate files.

### Required input files

#### Pipeline Configuration File (```-p```)
Modify example file to fit your needs. Primers can contain any base designation from the [IUPAC nucleotide code](https://www.bioinformatics.org/sms/iupac.html) in addition to I (inosine), which is considered as an "N" position. The recommendation for the ```blastLengthCutoff``` is 90% of the total marker target length. Suggested default ```taxonomyCutoffs``` are (97,95,90,80,70,60) for rRNA genes and (95,92,87,77,67,60) for protein encoding genes. The ```blastMode``` parameter has three options, with ```allEnvOUT``` and ```mostEnvOUT``` using negative taxID exclusion lists in the BLASTn step. ```mostEnvOUT``` removes subjects with "uncultured", "environmental samples", "metagenome", and "unidentified" in their "scientific name" field. ```allEnvOUT``` additionally removes subjects with "unclassified" in the "scientific name" field. Optional parameter ```failedMerge_useDirection``` should only be used if a prior run of the pipeline has shown that the forward/reverse reads will not merge properly. Optional parameter ```removeASVsFILE``` can be used to remove specified ASVs (one per line) from analysis (useful for known contaminants). Controls on DADA2 parameters are provided.

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
blastMode=mostEnvOUT #options: allIN,allEnvOUT,mostEnvOUT
```

#### Figure Configuration File (```-f```)
Modify example file to fit your needs. ```filterPercent``` controls the cutoff relative abundance that is used to place taxa into "zzOther" to simplify barchart figures. ```pieScale``` is an attempt to provide some control over the size of pie charts in the "mapPies"-labelled maps. ```removeNA``` removes any taxonomic assignment of "NA" from barchart figures. If ```filterLowQualSamples``` is "TRUE", then ```filterPercentLowQualSamples``` sets the cutoff for removing a sample from analysis, based on the average percentage of total reads of all other samples. Defining taxa of interest (```providedTaxaOfInterest``` "TRUE") can provide many additional figures focusing on only the listed taxa. The ```taxaOfInterestFile``` must contain a list of the taxa of interest (one per line), all from the same taxonomic level (in this case designated as "Order" for ```taxaOfInterestLevel```).

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

#### Sample Metadata File (```-s```)
This important file links a sample to its metadata. At a minimum, users must provide ```Sample``` designations, with ```lat``` and ```long``` recommended. ```Sample``` should be the first column; however, the order of other columns does not matter. Sample order (i.e. row order) is used for organizing figures. In addition to the controlled vocabulary (below), a user can provide an unlimited number of groupings (applied to visualizations, including NMDS/PCoA encircling) using the column headers from ```group1``` to ```group#```. Any columns that do not match the controlled vocabulary or ```group#``` are assigned to metadata and can be discrete or continuous variables (e.g. chemical/physical measurements, other observations). Any missing data should be filled in as ```NA```. Note that while phyloseq can handle missing data, vegan does not calculate evironmental fitting with missing data. We recommend that the pipeline is run first with all samples/data columns, followed up with a run to test environmental fitting with the samples with missing metadata and/or metadata columns with missing data removed.

Controlled vocabulary:
* ```Sample``` = Sample name identical to the ```X``` in ```X_R1.fastq.gz``` and ```X_R2.fastq.gz``` as found in the reads directory (```-r```). All sample names have ```MP_``` appended to the front and have non-alphanumeric characters converted to underscores to avoid programatic errors (e.g. R does not accept row/column headers that start with numbers).
* ```lat``` = latitude decimal degrees
* ```long``` = longitude decimal degrees
* ```replicates``` = String designating replicate samples. Replicate samples are treated both independently and as an averaged aggregate for some figures. Additional replicate-specific figures and tables are created, in particular for occupancy modelling outside the pipeline.
* ```sites``` = String designating samples taken from the same sites. Similar treatment to the ```replicates``` column in relation to averaged aggregate visualization.
* ```controls``` = Can be ```negative```, ```positive```, or ```NA```

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

#### Reads directory (```-r```)
Contains forward and reverse reads in the format ```X_R1.fastq.gz``` and ```X_R2.fastq.gz```, respectively, where ```X``` matches sample names in the sample metadata (```-s```) file. Fastq files must be gzipped.

### Optional input files

#### User-generated tab-delimited BLASTn file (```-b```)
After ASV-generation, a user can stop the pipeline to run their own BLASTn (use ```-outfmt '6 qseqid pident length staxids sacc'```) file on a custom database. To restart, delete the entries in ```progress.txt```, leaving only ```cutadaptFinished=TRUE``` and ```dada2_Finished=TRUE``` and restart the pipeline.

#### Import of SILVAngs results (```-e```)
This flag will trigger two prompts for file locations:
* "Enter the location of the SILVAngs ssu or lsu results directory (i.e. ~/Downloads/results/ssu)". This folder will contain the SILVAngs export file ```~/Downloads/results/ssu/exports/*---otus.csv``` and the SILVAngs cluster file ```~/Downloads/results/ssu/stats/sequence_cluster_map/data/*.fa.clstr```.
* "Enter the location of the reference taxonomy map for current SILVA database: i.e. tax_slv_ssu_138.1.txt". This file can be downloaded from the Arb SILVA [archive](https://www.arb-silva.de/no_cache/download/archive/current/Exports/taxonomy/).

## Stand alone applications

### morphology_convertTable2REVAMP.sh

### silvangs_convertTable2REVAMP.sh

### mergeBLASTandSILVAngs_ASV2Taxonomy.pl

### compare_markers_TablesFigures.sh

## REVAMP Results
TBD

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
