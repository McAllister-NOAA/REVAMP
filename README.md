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
* Krona Plots – [Ondov et al., 2011](https://doi.org/10.1186/1471-2105-12-385)
* phyloseq – [McMurdie and Holmes, 2013](https://doi.org/10.1371/journal.pone.0061217) ([useful tutorial](https://joey711.github.io/phyloseq/))
* vegan – [Oksanen et al., 2020](https://cran.r-project.org/web/packages/vegan/)

Many other auxiliary R packages make the figures possible: see dependencies.

## Wiki
REVAMP Wiki output guides, tutorials, and best practices is in development.

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

STOP HERE

NCBI databases (should be already installed; see above):
* *nt* database
* *taxonomy* database

Stand alone tools:
* R (Rscript) – v.4.0.3
* blast+ - v.2.13.0+
* cutadapt - v.2.8
* taxonkit - v.0.5.0
* TaxonomyTools (https://github.com/pmenzel/taxonomy-tools)
* KronaTools (https://github.com/marbl/Krona/wiki/Installing)
* subtree (taxonomy...)

R packages:
* dada2 – v.1.14.1
* dbplyr - v.1.4.2
* vegan - v.2.5-6
* mapping packages
* ggplot2 - v.3.3.0
* BiocManager
* phyloseq 
* ggpubr
* ggrepel
* ggalt
* mapdata
* mapproj
* marmap
* scatterpie
* janitor
* spatstat
* viridis
* raster
* rgdol

Perl dependencies:
* List::MoreUtils



## Toy Test dataset
TBD

## Quick Start Guide

### Required input files
* Raw reads
* Program settings
primer sequence F
primer sequence R
force merge? T/F [F]
DESeq rarefaction T/F [F]
expected insert size (bp)
sequence read length (bp)
location of sample metadata file
taxa depth (for figures)
control samples (positive/negative)
contaminant taxa list
taxid of interest list

* Sample metadata
sample order
sample groups
sample lat/long
replicate indication
chemistry

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
