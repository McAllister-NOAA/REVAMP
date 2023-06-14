#!/bin/bash

#NOTES: Require users to supply SILVAngs exports/x---[ls]su---otus.csv file from results archive.
#Pair with sample_metadata.txt for creating figures (same as metapipe metadata file).
## Sample Names should be identical to names in the sample metadata file.
## Sample Names will be cleaned of illegal characters (only alphanumeric and underline allowed).
## This cleaning step will also be applied to the sample names in the metadata file (so they are the same).
#
#NO spaces or quotes or other weird characters in file paths
#metapipe.sh must be in PATH

unset inputspreadsheet
unset silvaRefpath
unset clusterinfopath
unset samplemetafilepath
unset outdirectory
unset filterPercent
unset taxaOfInterestFile
unset taxaOfInterestCategory

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
sflag=0
oflag=0
fflag=0
rflag=0
bypassflag=FALSE
filterNAs=FALSE
providedTaxaOfInterest=FALSE
taxaOfInterestCategory=FALSE
mergeNCBISILVAeuks=FALSE
providedClusterInfo=FALSE

while getopts ":i:s:r:n:d:m:o:t:c:f:y" opt; do
  case ${opt} in
    i ) iflag=1
        inputspreadsheet=$OPTARG #Input SILVA spreadsheet
      ;;
    s ) sflag=1
        samplemetafilepath=$OPTARG #Sample metadata file
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
    m ) mergeNCBISILVAeuks=TRUE #Optional switch to create new folder with merged Eukaryotic assignments from NCBI with SILVA Bact/Arch assignments
      ;;
    r ) rflag=1
        silvaRefpath=$OPTARG #Input SILVA taxonomy map i.e. tax_slv_ssu_138.1.txt with tax\tsilvaAccession\tlevel info
      ;;
    d ) providedClusterInfo=TRUE
        clusterinfopath=$OPTARG #Optional cd-hit clustr mapping file
      ;;  
    y ) bypassflag=TRUE
      ;;
    \? ) echo "Invalid option: -$OPTARG"
         echo "Usage: silvangs_convertTable2Metapipe.sh" #Invalid option provided
         echo "       -i Input SILVAngs exports/x---[ls]su---otus.csv spreadsheet"
         echo "       -r Reference taxonomy map for current SILVA database: i.e. tax_slv_ssu_138.1.txt"
         echo "       -s Sample metadata file"
         echo "       -o Output directory"
         echo "       -f Filter percent cutoff for assignment to zzOther"
         echo "       -n Filter NAs from figures (optional)"
         echo "       -t Taxa of interest file (one per line) (optional)"
         echo "       -c Taxonomic category (e.g. Order) used in Taxa of interest file (required if -t called)"
         echo "       -m Merge NCBI Eukaryote taxa assignments with SILVA Bacteria/Archaea assignments (optional)"
         echo "       -d CD-HIT clustr mapping file if clustering done before hand (optional)"
         exit
      ;;
    : ) echo "Option is missing an argument: -$OPTARG"
        echo "Usage: silvangs_convertTable2Metapipe.sh" #Arg for a called option not provided
        echo "       -i Input SILVAngs exports/x---[ls]su---otus.csv spreadsheet"
        echo "       -r Reference taxonomy map for current SILVA database: i.e. tax_slv_ssu_138.1.txt"
        echo "       -s Sample metadata file"
        echo "       -o Output directory"
        echo "       -f Filter percent cutoff for assignment to zzOther"
        echo "       -n Filter NAs from figures (optional)"
        echo "       -t Taxa of interest file (one per line) (optional)"
        echo "       -c Taxonomic category (e.g. Order) used in Taxa of interest file (required if -t called)"
        echo "       -m Merge NCBI Eukaryote taxa assignments with SILVA Bacteria/Archaea assignments (optional)"
        echo "       -d CD-HIT clustr mapping file if clustering done before hand (optional)"
        exit
      ;;
  esac
done
shift $((OPTIND -1))

if [ $OPTIND -eq 1 ]
  then echo "Usage: silvangs_convertTable2Metapipe.sh" #No options passed
        echo "       -i Input SILVAngs exports/x---[ls]su---otus.csv spreadsheet"
        echo "       -r Reference taxonomy map for current SILVA database: i.e. tax_slv_ssu_138.1.txt"
        echo "       -s Sample metadata file"
        echo "       -o Output directory"
        echo "       -f Filter percent cutoff for assignment to zzOther"
        echo "       -n Filter NAs from figures (optional)"
        echo "       -t Taxa of interest file (one per line) (optional)"
        echo "       -c Taxonomic category (e.g. Order) used in Taxa of interest file (required if -t called)"
        echo "       -m Merge NCBI Eukaryote taxa assignments with SILVA Bacteria/Archaea assignments (optional)"
        echo "       -d CD-HIT clustr mapping file if clustering done before hand (optional)"
        exit
    fi

if [[ $iflag -eq 0 || $sflag -eq 0 || $oflag -eq 0 || $fflag -eq 0 || $rflag -eq 0 ]]
  then echo "All options except -n, -t, -c, -m, and -d are required."
        echo "Usage: silvangs_convertTable2Metapipe.sh" #Missing required options
        echo "       -i Input SILVAngs exports/x---[ls]su---otus.csv spreadsheet"
        echo "       -r Reference taxonomy map for current SILVA database: i.e. tax_slv_ssu_138.1.txt"
        echo "       -s Sample metadata file"
        echo "       -o Output directory"
        echo "       -f Filter percent cutoff for assignment to zzOther"
        echo "       -n Filter NAs from figures (optional)"
        echo "       -t Taxa of interest file (one per line) (optional)"
        echo "       -c Taxonomic category (e.g. Order) used in Taxa of interest file (required if -t called)"
        echo "       -m Merge NCBI Eukaryote taxa assignments with SILVA Bacteria/Archaea assignments (optional)"
        echo "       -d CD-HIT clustr mapping file if clustering done before hand (optional)"
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
perl ${metapipedir}/assets/sampleMetadata_fileCleanup.pl -i ${samplemetafilepath} > ${outdirectory}/sample_metadata_forR.txt

#Create ordered sample name file
cat ${outdirectory}/sample_metadata_forR.txt | cut -f1 | grep -v "Sample" > ${outdirectory}/sample_order.txt

##########################################################################################
##########################################################################################

##########################################################################################
##
##    Run perl script silvangs_convert.pl
##
##########################################################################################
if [[ "${providedClusterInfo}" = "TRUE" ]]; then
  if [[ "${mergeNCBISILVAeuks}" = "TRUE" ]]; then
    perl ${metapipedir}/assets/silvangs_convert.pl -i ${inputspreadsheet} -m ${metapipedir} -o ${outdirectory} -f ${filterPercent} -r ${silvaRefpath} -c ${clusterinfopath} -a
  else
    perl ${metapipedir}/assets/silvangs_convert.pl -i ${inputspreadsheet} -m ${metapipedir} -o ${outdirectory} -f ${filterPercent} -r ${silvaRefpath} -c ${clusterinfopath}
  fi
else
  if [[ "${mergeNCBISILVAeuks}" = "TRUE" ]]; then
    perl ${metapipedir}/assets/silvangs_convert.pl -i ${inputspreadsheet} -m ${metapipedir} -o ${outdirectory} -f ${filterPercent} -r ${silvaRefpath} -a
  else
    perl ${metapipedir}/assets/silvangs_convert.pl -i ${inputspreadsheet} -m ${metapipedir} -o ${outdirectory} -f ${filterPercent} -r ${silvaRefpath}
  fi
fi

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
while read -r line
  do mkdir -p ${outdirectory}/${line}_taxonomy/Figures/02_Barcharts/read_count
  mkdir -p ${outdirectory}/${line}_taxonomy/Figures/02_Barcharts/relative_abundance
  mkdir -p ${outdirectory}/${line}_taxonomy/Figures/03_Heatmaps/Taxonomy_merge_based
  mkdir -p ${outdirectory}/${line}_taxonomy/Figures/04_Alpha_Diversity/Taxonomy_merge_based
  mkdir -p ${outdirectory}/${line}_taxonomy/Figures/05_Ordination/Taxonomy_merge_based/read_count
  mkdir -p ${outdirectory}/${line}_taxonomy/Figures/05_Ordination/Taxonomy_merge_based/relative_abundance
  mkdir -p ${outdirectory}/${line}_taxonomy/Figures/05_Ordination/Taxonomy_merge_based/filterInclude_TOSPECIES_only/read_count
  mkdir -p ${outdirectory}/${line}_taxonomy/Figures/05_Ordination/Taxonomy_merge_based/filterInclude_TOSPECIES_only/relative_abundance
  mkdir -p ${outdirectory}/${line}_taxonomy/Figures/06_Network/Taxonomy_merge_based/read_count
  mkdir -p ${outdirectory}/${line}_taxonomy/Figures/06_Network/Taxonomy_merge_based/relative_abundance
  mkdir -p ${outdirectory}/${line}_taxonomy/Figures/07_Rarefaction_Curves
  mkdir -p ${outdirectory}/${line}_taxonomy/Figures/08_EnvironmentFit_Ordination/Taxonomy_merge_based
  
  if [[ "${providedTaxaOfInterest}" = "TRUE" ]]; then
  mkdir -p ${outdirectory}/${line}_taxonomy/Figures/Taxa_of_interest/02_Barcharts/read_count
  mkdir -p ${outdirectory}/${line}_taxonomy/Figures/Taxa_of_interest/02_Barcharts/relative_abundance
  mkdir -p ${outdirectory}/${line}_taxonomy/Figures/Taxa_of_interest/03_Heatmaps/Taxonomy_merge_based
  mkdir -p ${outdirectory}/${line}_taxonomy/Figures/Taxa_of_interest/06_Network/Taxonomy_merge_based/read_count
  mkdir -p ${outdirectory}/${line}_taxonomy/Figures/Taxa_of_interest/06_Network/Taxonomy_merge_based/relative_abundance
  fi
  
  #echo "${workingdirectory}/${outdirectory}/${line}_taxonomy/Figures ${workingdirectory}/${outdirectory}/${line}_taxonomy/silvangs_${line}_asvTaxonomyTable_NOUNKNOWNS.txt ${workingdirectory}/${outdirectory}/${line}_taxonomy/ASVs_counts_NOUNKNOWNS.tsv ${workingdirectory}/${outdirectory}/${line}_taxonomy/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.tsv ${outdirectory} ${workingdirectory}/${outdirectory}/sample_metadata_forR.txt $replicates $sites $filterNAs $groupsDefinedFlag $numberGroupsDefined $filterNAs $providedTaxaOfInterest $taxaOfInterestCategory ${workingdirectory}/${outdirectory}/taxaOfInterest.txt $chemData $locationChemHeaders ${workingdirectory}/${outdirectory}/sample_order.txt ${workingdirectory}/${outdirectory}/${line}_taxonomy/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt $filterPercent"
  Rscript --vanilla ${metapipedir}/assets/phyloseq_collapsedOnTaxonomy_individualRun.R ${workingdirectory}/${outdirectory}/${line}_taxonomy/Figures ${workingdirectory}/${outdirectory}/${line}_taxonomy/silvangs_${line}_asvTaxonomyTable_NOUNKNOWNS.txt ${workingdirectory}/${outdirectory}/${line}_taxonomy/ASVs_counts_NOUNKNOWNS.tsv ${workingdirectory}/${outdirectory}/${line}_taxonomy/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.tsv ${outdirectory} ${workingdirectory}/${outdirectory}/sample_metadata_forR.txt $replicates $sites $filterNAs $groupsDefinedFlag $numberGroupsDefined $filterNAs $providedTaxaOfInterest $taxaOfInterestCategory ${workingdirectory}/${outdirectory}/taxaOfInterest.txt $chemData $locationChemHeaders ${workingdirectory}/${outdirectory}/sample_order.txt ${workingdirectory}/${outdirectory}/${line}_taxonomy/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt $filterPercent \
    1>> ${workingdirectory}/${outdirectory}/${line}_taxonomy/Figures/phyloseq_rscript_out.log 2>&1
  #rm -f ${workingdirectory}/${outdirectory}/Figures/02_Barcharts/read_count/Rplots.pdf
  
  Rscript --vanilla ${metapipedir}/assets/environment_fit_ordination_collapsedOnTaxonomy_individualrun.R ${workingdirectory}/${outdirectory}/${line}_taxonomy/Figures ${workingdirectory}/${outdirectory}/${line}_taxonomy/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.tsv ${workingdirectory}/${outdirectory}/sample_metadata_forR.txt $replicates $sites $chemData $locationChemHeaders \
    1>> ${workingdirectory}/${outdirectory}/${line}_taxonomy/Figures/envfit_rscript_out.log 2>&1
  
done < ${outdirectory}/unique_TaxaFolders.txt

echo "YOU MADE IT!"