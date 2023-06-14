#!/bin/bash

#NOTES:Requires FILE1) ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.txt, FILE2) x_asvTaxonomyTable_NOUNKNOWNS.txt, FILE3) sample_metadata_forR.txt, FILE4) taxonomy2PercentAbundance_humanReadable.txt
# Copy these files into a single folder (-i) with the following four folders within it:
# "ASV_relabund" "ASV_taxonomy" "Sample_metadata" "Taxa_relabund_Human"
# Rename the files in each folder to be "Marker.txt" (i.e. "COI.txt").
# Will also need a file for sample equivalents with no header: "Marker(no .txt)\tMarkerSample\tReferenceSample" (human readable comparisons to reference marker/method only).
# Note: If the "group#" do not align between marker sample metadata files, those groups will not match across markers and resulting figures based on groups should be ignored.

#NO spaces or quotes or other weird characters in file paths
#metapipe.sh must be in PATH

unset inputfolder
unset outdirectory
unset taxaOfInterestFile
unset taxaOfInterestCategory
unset sampleEquivalents
unset refMarker

workingdirectory=`pwd`

unset metapipedir
unset tempprogdir
tempprogdir=`which metapipe.sh`
metapipedir=`echo $tempprogdir | sed -E 's/\/metapipe\.sh$//'`
myInvocation="$(printf %q "$BASH_SOURCE")$((($#)) && printf ' %q' "$@")"

##########################################################################################
##
##    Get command line options
##
##########################################################################################
iflag=0
oflag=0
providedTaxaOfInterest=FALSE
taxaOfInterestCategory=FALSE
runHumanReadableComparisons=FALSE

while getopts ":i:o:s:r:t:c:" opt; do
  case ${opt} in
    i ) iflag=1
        inputdirectory=$OPTARG #Input folder with four folders as described above
        inputdirectory=`echo $inputdirectory | sed -E 's/\/$//'`
      ;;
    o ) oflag=1
        outdirectory=$OPTARG #Location for output files
        outdirectory=`echo $outdirectory | sed -E 's/\/$//'`
      ;;
    s ) runHumanReadableComparisons=TRUE
        sampleEquivalents=$OPTARG #Location of sample equivalents file
      ;;
    r ) refMarker=$OPTARG #Choice of marker as reference in comparison (marker.txt)
      ;;
    t ) providedTaxaOfInterest=TRUE
        taxaOfInterestFile=$OPTARG #Location of taxa of interest file
      ;;
    c ) taxaOfInterestCategory=$OPTARG #Taxonomic category represented in the Taxa of Interest File
      ;;
    \? ) echo "Invalid option: -$OPTARG"
         echo "Usage: compare_markers_TablesFigures.sh" #Invalid option provided
         echo "       -i Input folder with internal folders: ASV_relabund (req), ASV_taxonomy (req), Sample_metadata (req), and Taxa_relabund_Human (opt)"
         echo "          Each folder should have their corresponding files renamed to Marker.txt (i.e. COI.txt)"
         echo "       -o Output directory"
         echo "       -s Sample equivalents file (Marker tab MarkerSample tab ReferenceSample) (optional)"
         echo "       -r Reference marker file name (Marker.txt) (required if -s is called)"
         echo "       -t Taxa of interest file (one per line) (optional)"
         echo "       -c Taxonomic category (e.g. Order) used in Taxa of interest file (required if -t called)"
         exit
      ;;
    : ) echo "Option is missing an argument: -$OPTARG"
        echo "Usage: compare_markers_TablesFigures.sh" #Arg for a called option not provided
        echo "       -i Input folder with internal folders: ASV_relabund (req), ASV_taxonomy (req), Sample_metadata (req), and Taxa_relabund_Human (opt)"
        echo "          Each folder should have their corresponding files renamed to Marker.txt (i.e. COI.txt)"
        echo "       -o Output directory"
        echo "       -s Sample equivalents file (Marker tab MarkerSample tab ReferenceSample) (optional)"
        echo "       -r Reference marker file name (Marker.txt) (required if -s is called)"
        echo "       -t Taxa of interest file (one per line) (optional)"
        echo "       -c Taxonomic category (e.g. Order) used in Taxa of interest file (required if -t called)"
        exit
      ;;
  esac
done
shift $((OPTIND -1))

if [ $OPTIND -eq 1 ]
  then echo "Usage: compare_markers_TablesFigures.sh" #No options passed
       echo "       -i Input folder with internal folders: ASV_relabund (req), ASV_taxonomy (req), Sample_metadata (req), and Taxa_relabund_Human (opt)"
       echo "          Each folder should have their corresponding files renamed to Marker.txt (i.e. COI.txt)"
       echo "       -o Output directory"
       echo "       -s Sample equivalents file (Marker tab MarkerSample tab ReferenceSample) (optional)"
       echo "       -r Reference marker file name (Marker.txt) (required if -s is called)"
       echo "       -t Taxa of interest file (one per line) (optional)"
       echo "       -c Taxonomic category (e.g. Order) used in Taxa of interest file (required if -t called)"
       exit
    fi

if [[ $iflag -eq 0 || $oflag -eq 0 ]]
  then echo "All options except -s, -r, -t, and -c are required."
        echo "Usage: compare_markers_TablesFigures.sh" #Missing required options
        echo "       -i Input folder with internal folders: ASV_relabund (req), ASV_taxonomy (req), Sample_metadata (req), and Taxa_relabund_Human (opt)"
        echo "          Each folder should have their corresponding files renamed to Marker.txt (i.e. COI.txt)"
        echo "       -o Output directory"
        echo "       -s Sample equivalents file (Marker tab MarkerSample tab ReferenceSample) (optional)"
        echo "       -r Reference marker file name (Marker.txt) (required if -s is called)"
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
mkdir -p ${outdirectory}/Tables
mkdir -p ${outdirectory}/Figures/Ordination
mkdir -p ${outdirectory}/Figures/Network
mkdir -p ${outdirectory}/Figures/VennDiagrams


cp -r $inputdirectory ${outdirectory}/comparisons_inputDir

touch ${outdirectory}/run.log

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

#Make area for human readable comparisons
if [[ "${runHumanReadableComparisons}" = "TRUE" ]]; then
  mkdir -p ${outdirectory}/Tables/HumanReadable
  cat ${sampleEquivalents} | sed "s/$(printf '\r')\$//" > ${outdirectory}/Tables/HumanReadable/sampleEquivalents.txt
fi

##########################################################################################
##########################################################################################

##########################################################################################
##
##    Run perl script compare_markers.pl
##
##########################################################################################
if [[ "${providedTaxaOfInterest}" = "TRUE" ]]; then
  perl ${metapipedir}/assets/compare_markers.pl -i ${inputdirectory}/ASV_taxonomy -o ${outdirectory}/Tables -f ${outdirectory}/taxaOfInterest.txt -l ${taxaOfInterestCategory}
else
  perl ${metapipedir}/assets/compare_markers.pl -i ${inputdirectory}/ASV_taxonomy -o ${outdirectory}/Tables
fi

##########################################################################################
##########################################################################################

##########################################################################################
##
##    Run perl script compare_marker_cleanup.pl in preparation for R phyloseq script
##
##########################################################################################
perl ${metapipedir}/assets/compare_marker_cleanup.pl -i $inputdirectory -o $outdirectory

##########################################################################################
##########################################################################################

##########################################################################################
##
##    Run perl script compare_humanReadableTables.pl to create comparison table (optional)
##
##########################################################################################
if [[ "${runHumanReadableComparisons}" = "TRUE" ]]; then
  perl ${metapipedir}/assets/compare_humanReadableTables.pl -i ${inputdirectory}/Taxa_relabund_Human -r $refMarker -s ${outdirectory}/Tables/HumanReadable/sampleEquivalents.txt > ${outdirectory}/Tables/HumanReadable/comparison_HumanReadableTables_2_Reference.txt
fi

##########################################################################################
##########################################################################################

##########################################################################################
##
##    FIGURES
##
##########################################################################################
#Interpreted from sample metadata file:
replicates="FALSE"
sites="FALSE"
groupsDefinedFlag="FALSE"

#Interpret sample metadata file for inputs:
highestgroupnum=0

for ((f=1; f<=`awk '{print NF}' ${workingdirectory}/${outdirectory}/MergedMarkers_sample_metadata_forR.txt | sort -nu | tail -n 1`; f++))
  do cat ${workingdirectory}/${outdirectory}/MergedMarkers_sample_metadata_forR.txt | cut -f${f} > ${workingdirectory}/${outdirectory}/temp
    header=`head -n 1 ${workingdirectory}/${outdirectory}/temp`
    if [[ "$header" = "Sample" ]]; then
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
    fi
  done
  rm ${workingdirectory}/${outdirectory}/temp
  numberGroupsDefined=$highestgroupnum

#Phyloseq figures
Rscript --vanilla ${metapipedir}/assets/phyloseq_mergedMarkerComparison.R ${workingdirectory}/${outdirectory} ${workingdirectory}/${outdirectory}/MergedMarkers_asvTaxonomyTable_NOUNKNOWNS.txt ${workingdirectory}/${outdirectory}/MergedMarkers_ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.txt ${workingdirectory}/${outdirectory}/MergedMarkers_sample_metadata_forR.txt $replicates $sites $groupsDefinedFlag $numberGroupsDefined $providedTaxaOfInterest $taxaOfInterestCategory ${workingdirectory}/${outdirectory}/taxaOfInterest.txt \
    1>> ${workingdirectory}/${outdirectory}/phyloseq_rscript_out.log 2>&1

rm -f ${workingdirectory}/${outdirectory}/Figures/Ordination/Rplots.pdf

#Venn diagrams (only if there are 2-5 marker comparisons)
numVennComp=`ls ${inputdirectory}/ASV_relabund | wc -l`

declare -a hierarchyArr=("Kingdom" "Phylum" "Class" "Order" "Family" "Genus" "Species")
declare -a uniqueTotTerm=("Total" "Terminal")

for i in "${uniqueTotTerm[@]}"
do
  for j in "${hierarchyArr[@]}"
  do
    Rscript --vanilla ${metapipedir}/assets/venn_diagrams.R ${workingdirectory}/${outdirectory} $numVennComp $j $i \
      1>> ${workingdirectory}/${outdirectory}/venn_rscript_out.log 2>&1
  done
done

rm -f ${workingdirectory}/${outdirectory}/Figures/VennDiagrams/Rplots.pdf

if [[ $numVennComp -lt 2 || $numVennComp -gt 5 ]]; then
  rmdir ${workingdirectory}/${outdirectory}/Figures/VennDiagrams
fi

echo "YOU MADE IT!"
