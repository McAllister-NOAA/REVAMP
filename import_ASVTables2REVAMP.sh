#!/bin/bash

#NOTES: Require users to supply base ASV count table, taxa table, percent abundance table and sample metadata table.
#Eventually deal with zzOther TBD
## Names should be identical to names in the sample metadata file.
## Names will NOT be cleaned of illegal characters (only alphanumeric and underline allowed). So sample metadata and ASV counts should already match.
#
#NO spaces or quotes or other weird characters in file paths
#revamp.sh must be in PATH
#Will need taxonomydmp - run prepscript beforehand

unset inputspreadsheet
unset taxonomyfilepath
unset percabundpath
unset samplemetafilepath
unset outdirectory
unset filterPercent
unset taxaOfInterestFile
unset taxaOfInterestCategory

workingdirectory=`pwd`

unset revampdir
unset tempprogdir
tempprogdir=`which revamp.sh`
revampdir=`echo $tempprogdir | sed -E 's/\/revamp\.sh$//'`
myInvocation="$(printf %q "$BASH_SOURCE")$((($#)) && printf ' %q' "$@")"

##########################################################################################
##
##    Get command line options
##
##########################################################################################
iflag=0
sflag=0
oflag=0
fflag=0
xflag=0
pflag=0
bypassflag=FALSE
filterNAs=FALSE
providedTaxaOfInterest=FALSE
taxaOfInterestCategory=FALSE

while getopts ":i:s:x:p:n:o:t:c:f:y" opt; do
  case ${opt} in
    i ) iflag=1
        inputspreadsheet=$OPTARG #Input ASV count file
      ;;
    s ) sflag=1
        samplemetafilepath=$OPTARG #Sample metadata file
      ;;
    x ) xflag=1
        taxonomyfilepath=$OPTARG #ASV Taxonomy file
      ;;
    p ) pflag=1
        percabundpath=$OPTARG #ASV percent abundance file
      ;;
    n ) filterNAs=TRUE #Whether to remove NAs from phyloseq figures
      ;;
    o ) oflag=1
        outdirectory=$OPTARG #Location for output files
        outdirectory=`echo $outdirectory | sed -E 's/\/$//'`
      ;;
    f ) fflag=1
        filterPercent=$OPTARG #Percent filter for taxa to be assigned to zzOther
      ;;
    t ) providedTaxaOfInterest=TRUE
        taxaOfInterestFile=$OPTARG #Location of taxa of interest file
      ;;
    c ) taxaOfInterestCategory=$OPTARG #Taxonomic category represented in the Taxa of Interest File
      ;;
    y ) bypassflag=TRUE
      ;;
    \? ) echo "Invalid option: -$OPTARG"
         echo "Usage: morphology_convertTable2REVAMP.sh" #Invalid option provided
         echo "       -i ASV count file"
         echo "       -s Sample metadata file"
         echo "       -x ASV Taxonomy file"
         echo "       -p ASV percent abundance table"
         echo "       -o Output directory"
         echo "       -f Filter percent cutoff for assignment to zzOther"
         echo "       -n Filter NAs from figures (optional)"
         echo "       -t Taxa of interest file (one per line) (optional)"
         echo "       -c Taxonomic category (e.g. Order) used in Taxa of interest file (required if -t called)"
         exit
      ;;
    : ) echo "Option is missing an argument: -$OPTARG"
        echo "Usage: morphology_convertTable2REVAMP.sh" #Arg for a called option not provided
        echo "       -i ASV count file"
        echo "       -s Sample metadata file"
        echo "       -x ASV Taxonomy file"
        echo "       -p ASV percent abundance table"
        echo "       -o Output directory"
        echo "       -f Filter percent cutoff for assignment to zzOther"
        echo "       -n Filter NAs from figures (optional)"
        echo "       -t Taxa of interest file (one per line) (optional)"
        echo "       -c Taxonomic category (e.g. Order) used in Taxa of interest file (required if -t called)"
        exit
      ;;
  esac
done
shift $((OPTIND -1))

if [ $OPTIND -eq 1 ]
  then echo "Usage: morphology_convertTable2REVAMP.sh" #No options passed
        echo "       -i ASV count file"
        echo "       -s Sample metadata file"
        echo "       -x ASV Taxonomy file"
        echo "       -p ASV percent abundance table"
        echo "       -o Output directory"
        echo "       -f Filter percent cutoff for assignment to zzOther"
        echo "       -n Filter NAs from figures (optional)"
        echo "       -t Taxa of interest file (one per line) (optional)"
        echo "       -c Taxonomic category (e.g. Order) used in Taxa of interest file (required if -t called)"
        exit
    fi

if [[ $iflag -eq 0 || $pflag -eq 0 || $xflag -eq 0 || $sflag -eq 0 || $oflag -eq 0 || $fflag -eq 0 ]]
  then echo "All options except -n, -t, and -c are required."
        echo "Usage: morphology_convertTable2REVAMP.sh" #Missing required options
        echo "       -i ASV count file"
        echo "       -s Sample metadata file"
        echo "       -x ASV Taxonomy file"
        echo "       -p ASV percent abundance table"
        echo "       -o Output directory"
        echo "       -f Filter percent cutoff for assignment to zzOther"
        echo "       -n Filter NAs from figures (optional)"
        echo "       -t Taxa of interest file (one per line) (optional)"
        echo "       -c Taxonomic category (e.g. Order) used in Taxa of interest file (required if -t called)"
        exit
      fi
##########################################################################################
##########################################################################################

##########################################################################################
##
##    Cleanup and prep
##
##########################################################################################
mkdir ${outdirectory}
touch ${outdirectory}/run.log
cp ${samplemetafilepath} ${outdirectory}/sample_metadata.txt
cp ${inputspreadsheet} ${outdirectory}/ASVs_counts_used.tsv
cp ${taxonomyfilepath} ${outdirectory}/ASVTaxonomyTable_used.txt
cp ${percabundpath} ${outdirectory}/ASVs_percAbundance_used.tsv

exec &> >(tee -a ${outdirectory}/run.log)

echo
echo "Start of run:"
date
echo
echo "Invoked script options:"
echo "$myInvocation"
echo

#Copy taxa of interest file to outdirectory:
if [[ "${providedTaxaOfInterest}" = "TRUE" ]]; then
  cp $taxaOfInterestFile ${outdirectory}/taxaOfInterest.txt
  fi

#Create sample metadata file with identical manipulation of sample names for downstream R work
perl ${revampdir}/assets/sampleMetadata_fileCleanup.pl -i ${samplemetafilepath} > ${outdirectory}/sample_metadata_forR.txt

cat ${outdirectory}/ASVs_counts_used.tsv | sed -E '1 s/[^a-zA-Z0-9_[:space:]]/_/g' > ${outdirectory}/ASVs_counts_used_forR.tsv
cat ${outdirectory}/ASVs_percAbundance_used.tsv | sed -E '1 s/[^a-zA-Z0-9_[:space:]]/_/g' > ${outdirectory}/ASVs_percAbundance_used_forR.tsv

#Create ordered sample name file
cat ${outdirectory}/sample_metadata_forR.txt | cut -f1 | grep -v "Sample" > ${outdirectory}/sample_order.txt

##########################################################################################
##########################################################################################

##########################################################################################
##
##    FIGURES & ANALYSIS
##
##########################################################################################
#Interpreted from sample metadata file:
replicates="FALSE"
sites="FALSE"
chemData="FALSE"
groupsDefinedFlag="FALSE"
locationChemHeaders="NULL"


#Interpret sample metadata file for inputs:
highestgroupnum=0
rm -f ${workingdirectory}/${outdirectory}/chem_headers.txt

for ((f=1; f<=`awk '{print NF}' ${workingdirectory}/${outdirectory}/sample_metadata_forR.txt | sort -nu | tail -n 1`; f++))
  do cat ${workingdirectory}/${outdirectory}/sample_metadata_forR.txt | cut -f${f} > ${workingdirectory}/${outdirectory}/temp
    header=`head -n 1 ${workingdirectory}/${outdirectory}/temp`
    if [[ "$header" = "Sample" || "$header" = "lat" || "$header" = "long" ]]; then
      continue
      elif [[ "$header" = "replicates" ]]; then
      replicates="TRUE"
      elif [[ "$header" = "sites" ]]; then
      sites="TRUE"
      elif [[ "$header" =~ "group" ]]; then
      groupsDefinedFlag="TRUE"
      numGrp=`echo $header | sed -E 's/group//'`
      if [[ $numGrp -gt $highestgroupnum ]]; then
        highestgroupnum=$numGrp
      fi
    else
      chemData="TRUE"
      echo $header >> ${workingdirectory}/${outdirectory}/chem_headers.txt
      locationChemHeaders="${workingdirectory}/${outdirectory}/chem_headers.txt"
    fi
  done
  rm ${workingdirectory}/${outdirectory}/temp
  numberGroupsDefined=$highestgroupnum

#Phyloseq figures
mkdir -p ${outdirectory}/Figures/02_Barcharts/read_count
mkdir -p ${outdirectory}/Figures/02_Barcharts/relative_abundance
mkdir -p ${outdirectory}/Figures/03_Heatmaps/Taxonomy_merge_based
mkdir -p ${outdirectory}/Figures/04_Alpha_Diversity/Taxonomy_merge_based
mkdir -p ${outdirectory}/Figures/05_Ordination/Taxonomy_merge_based/read_count
mkdir -p ${outdirectory}/Figures/05_Ordination/Taxonomy_merge_based/relative_abundance
mkdir -p ${outdirectory}/Figures/05_Ordination/Taxonomy_merge_based/filterInclude_TOSPECIES_only/read_count
mkdir -p ${outdirectory}/Figures/05_Ordination/Taxonomy_merge_based/filterInclude_TOSPECIES_only/relative_abundance
mkdir -p ${outdirectory}/Figures/06_Network/Taxonomy_merge_based/read_count
mkdir -p ${outdirectory}/Figures/06_Network/Taxonomy_merge_based/relative_abundance
mkdir -p ${outdirectory}/Figures/07_Rarefaction_Curves
mkdir -p ${outdirectory}/Figures/08_EnvironmentFit_Ordination/Taxonomy_merge_based

if [[ "${providedTaxaOfInterest}" = "TRUE" ]]; then
mkdir -p ${outdirectory}/Figures/Taxa_of_interest/02_Barcharts/read_count
mkdir -p ${outdirectory}/Figures/Taxa_of_interest/02_Barcharts/relative_abundance
mkdir -p ${outdirectory}/Figures/Taxa_of_interest/03_Heatmaps/Taxonomy_merge_based
mkdir -p ${outdirectory}/Figures/Taxa_of_interest/06_Network/Taxonomy_merge_based/read_count
mkdir -p ${outdirectory}/Figures/Taxa_of_interest/06_Network/Taxonomy_merge_based/relative_abundance
fi

Rscript --vanilla ${revampdir}/assets/phyloseq_collapsedOnTaxonomy_individualRun.R ${workingdirectory}/${outdirectory}/Figures ${workingdirectory}/${outdirectory}/ASVTaxonomyTable_used.txt ${workingdirectory}/${outdirectory}/ASVs_counts_used_forR.tsv ${workingdirectory}/${outdirectory}/ASVs_percAbundance_used_forR.tsv ${outdirectory} ${workingdirectory}/${outdirectory}/sample_metadata_forR.txt $replicates $sites $filterNAs $groupsDefinedFlag $numberGroupsDefined $filterNAs $providedTaxaOfInterest $taxaOfInterestCategory ${workingdirectory}/${outdirectory}/taxaOfInterest.txt $chemData $locationChemHeaders ${workingdirectory}/${outdirectory}/sample_order.txt ${workingdirectory}/${outdirectory}/ASVTaxonomyTable_used.txt $filterPercent \
  1>> ${workingdirectory}/${outdirectory}/Figures/phyloseq_rscript_out.log 2>&1
#rm -f ${workingdirectory}/${outdirectory}/Figures/02_Barcharts/read_count/Rplots.pdf

Rscript --vanilla ${revampdir}/assets/environment_fit_ordination_collapsedOnTaxonomy_individualrun.R ${workingdirectory}/${outdirectory}/Figures ${workingdirectory}/${outdirectory}/ASVs_counts_used_forR.tsv ${workingdirectory}/${outdirectory}/sample_metadata_forR.txt $replicates $sites $chemData $locationChemHeaders \
  1>> ${workingdirectory}/${outdirectory}/Figures/envfit_rscript_out.log 2>&1
  


echo "YOU MADE IT!"
