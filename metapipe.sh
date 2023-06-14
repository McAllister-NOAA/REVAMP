#!/bin/bash

#NOTES: Require users to supply raw sequence data files in format "samplename_R[12].fastq.gz".
## These names should be identical to names in the sample metadata file.
## Names will be cleaned of illegal characters (only alphanumeric and underline allowed).
## This cleaning step will also be applied to the sample names in the metadata file (so they are the same).
#
#NO spaces or quotes or other weird characters in file paths
#metapipe.sh must be in PATH
#Pro tip: steps that use memory only use 70% of what max you give
#Will need nt and taxonomydmp - run prepscript beforehand
#Can provide own blastn output on whatever database. Must be formatted -outfmt '6 qseqid pident length staxids'

unset parameterfilepath
unset samplemetafilepath
unset readfolderpath
unset outdirectory
unset optionalUserBLASTResult
unset figureparamfilepath
unset silvangsInputExportFile
unset silvangsInputClusterFile
unset silvaRefMap

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
pflag=0
sflag=0
rflag=0
oflag=0
fflag=0
cflag=0
blastflag=FALSE
silvaASVflag=FALSE
bypassflag=FALSE
keepIntermediateFiles=TRUE

while getopts ":p:s:r:o:b:f:yke" opt; do
  case ${opt} in
    p ) pflag=1
        parameterfilepath=$OPTARG #metapipe_config.txt (see README)
      ;;
    s ) sflag=1
        samplemetafilepath=$OPTARG #Sample metadata file (see README)
      ;;
    r ) rflag=1
        readfolderpath=$OPTARG #Folder with fastq.gz raw read files
      ;;
    o ) oflag=1
        outdirectory=$OPTARG #Location for output files
        outdirectory=`echo $outdirectory | sed -E 's/\/$//'`
      ;;
    f ) fflag=1
        figureparamfilepath=$OPTARG #Location of figure config file (see README)
      ;;
    b ) blastflag=TRUE
        optionalUserBLASTResult=$OPTARG #Location of user blastn input (optional)
      ;;
    e ) silvaASVflag=TRUE
      ;;
    y ) bypassflag=TRUE
      ;;
    k ) keepIntermediateFiles=FALSE
      ;;
    \? ) echo "Invalid option: -$OPTARG"
         echo "Usage: metapipe.sh" #Invalid option provided
         echo "       -p Config File"
         echo "       -f Figure config file"
         echo "       -s Sample metadata file"
         echo "       -r Read folder"
         echo "       -o Output directory"
         echo "       -b User-supplied BLASTn btab result file (optional)"
         echo "       -e Toggle use of SILVAngs taxonomy assignments by ASV (optional)"
         echo "       -y Bypass all terminal prompts (optional)"
         echo "       -k Remove all intermediate files (optional)"
         echo "          Not recommended for first run"
         echo "          (best to rerun with -k after successful completion)"
         echo "       See README for details"
         exit
      ;;
    : ) echo "Option is missing an argument: -$OPTARG"
        echo "Usage: metapipe.sh" #Arg for a called option not provided
        echo "       -p Config File"
        echo "       -f Figure config file"
        echo "       -s Sample metadata file"
        echo "       -r Read folder"
        echo "       -o Output directory"
        echo "       -b User-supplied BLASTn btab result file (optional)"
        echo "       -e Toggle use of SILVAngs taxonomy assignments by ASV (optional)"
        echo "       -y Bypass all terminal prompts (optional)"
        echo "       -k Remove all intermediate files (optional)"
        echo "          Not recommended for first run"
        echo "          (best to rerun with -k after successful completion)"
        echo "       See README for details"
        exit
      ;;
  esac
done
shift $((OPTIND -1))

if [ $OPTIND -eq 1 ]
  then echo "Usage: metapipe.sh" #No options passed
        echo "       -p Config File"
        echo "       -f Figure config file"
        echo "       -s Sample metadata file"
        echo "       -r Read folder"
        echo "       -o Output directory"
        echo "       -b User-supplied BLASTn btab result file (optional)"
        echo "       -e Toggle use of SILVAngs taxonomy assignments by ASV (optional)"
        echo "       -y Bypass all terminal prompts (optional)"
        echo "       -k Remove all intermediate files (optional)"
        echo "          Not recommended for first run"
        echo "          (best to rerun with -k after successful completion)"
        echo "       See README for details"
        exit
    fi

if [[ $pflag -eq 0 || $sflag -eq 0 || $rflag -eq 0 || $oflag -eq 0 || $fflag -eq 0 ]]
  then echo "All options except -b, -y, and -k are required."
        echo "Usage: metapipe.sh" #Missing required options
        echo "       -p Config File"
        echo "       -f Figure config file"
        echo "       -s Sample metadata file"
        echo "       -r Read folder"
        echo "       -o Output directory"
        echo "       -b User-supplied BLASTn btab result file (optional)"
        echo "       -e Toggle use of SILVAngs taxonomy assignments by ASV (optional)"
        echo "       -y Bypass all terminal prompts (optional)"
        echo "       -k Remove all intermediate files (optional)"
        echo "          Not recommended for first run"
        echo "          (best to rerun with -k after successful completion)"
        echo "       See README for details"
        exit
      fi
      
##########################################################################################
##########################################################################################

##########################################################################################
##
##    Get program config file options
##
##########################################################################################
unset primerF
unset primerR
unset blastLengthCutoff
unset systemmemoryMB
unset locationNTdatabase
unset taxonomyCutoffs
unset dada_minlength
unset dada_phix
unset dada_trunQ
unset dada_maxEE1
unset dada_maxEE2
unset dada_trimRight
unset dada_trimLeft
unset blastMode
failedMerge_useDirection=NA
removeASVsFILE=NULL

source $parameterfilepath

temp_holder=$taxonomyCutoffs
taxonomyCutoffs=`echo $temp_holder | sed -E 's/[^0-9]/,/g'`

#TEMP Holding:
#echo "<<<$taxonomyCutoffs>>>"
#sed -E 's/^\D//' | sed -E 's/\D$//' 

##########################################################################################
##########################################################################################

##########################################################################################
##
##    Cleanup and prep
##
##########################################################################################
cutadaptFinished=FALSE
dada2_Finished=FALSE
blastFinished=FALSE
blastformattingFinished=FALSE
taxonomyscriptFinished=FALSE
figuresFinished=FALSE

if [ -d "${outdirectory}" ]; then
  configcompare1=`cat ${parameterfilepath}`
  configcompare2=`cat ${outdirectory}/config_file.txt`
  samplecompare1=`cat ${samplemetafilepath}`
  samplecompare2=`cat ${outdirectory}/sample_metadata.txt`
  if [[ "${configcompare1}" = "${configcompare2}" ]]; then
    echo "Config files identical"
  else
    echo "Config files differ between runs, choose a different out directory"
    exit
  fi
  if [[ "${samplecompare1}" = "${samplecompare2}" ]]; then
    echo "Sample metadata files identical"
  else
    echo "Sample metadata files differ between runs, choose a different out directory"
    exit
  fi
  source ${outdirectory}/progress.txt
  mkdir -p ${outdirectory}/run_logs
  currenttime=`date | sed -E 's/[^A-Za-z0-9_]/_/g'`
  mv ${outdirectory}/run.log ${outdirectory}/run_logs/runlog_${currenttime}.txt
  touch ${outdirectory}/run.log
  echo >> ${outdirectory}/Rscript_arguments.log
  echo "New Run" >> ${outdirectory}/Rscript_arguments.log
else
  mkdir ${outdirectory}
  touch ${outdirectory}/progress.txt
  touch ${outdirectory}/run.log
  echo "First Run" >> ${outdirectory}/Rscript_arguments.log
  cp ${parameterfilepath} ${outdirectory}/config_file.txt
  cp ${samplemetafilepath} ${outdirectory}/sample_metadata.txt
fi

exec &> >(tee -a ${outdirectory}/run.log)

echo
echo "Start of run:"
date
echo
echo "Invoked script options:"
echo "$myInvocation"
echo

#Clean raw read files to remove illegal characters
cd ${readfolderpath}
for f in *fastq.gz
do cleanprior=`echo $f | sed -E 's/^MP_//'`
        base=$(basename $cleanprior .fastq.gz)
        newname=`echo ${base} | sed -E 's/[^A-Za-z0-9_]/_/g'`
        mv $f MP_${newname}.fastq.gz
        done
cd ${workingdirectory}

#Create sample metadata file with identical manipulation of sample names for downstream R work
perl ${metapipedir}/assets/sampleMetadata_fileCleanup.pl -i ${samplemetafilepath} > ${outdirectory}/sample_metadata_forR.txt

#Create ordered sample name file
cat ${outdirectory}/sample_metadata_forR.txt | cut -f1 | grep -v "Sample" > ${outdirectory}/sample_order.txt


##########################################################################################
##########################################################################################

##########################################################################################
##
##    Run CUTADAPT (v.3.7)
##
##########################################################################################
if [[ "${cutadaptFinished}" = "TRUE" ]]; then
  echo "Cutadapt from prior run"
else
  echo "Running Cutadapt: `date`"
  if [ -d "${outdirectory}/cutadapt" ]; then
    rm -r ${outdirectory}/cutadapt
    mkdir ${outdirectory}/cutadapt
  else
    mkdir ${outdirectory}/cutadapt
  fi
    
  #Reverse complement primers (I is converted to N; T/U ok)
  revcomp_primerF=`echo $primerF | tr ACGTUWSMKRYBDHVNIacgtuwsmkrybdhvni TGCAAWSKMYRVHDBNNtgcaawskmyrvhdbnn | rev`
  revcomp_primerR=`echo $primerR | tr ACGTUWSMKRYBDHVNIacgtuwsmkrybdhvni TGCAAWSKMYRVHDBNNtgcaawskmyrvhdbnn | rev`
  
  temp_primerF=`echo $primerF | tr Ii Nn`
  temp_primerR=`echo $primerR | tr Ii Nn`
  primerF=$temp_primerF
  primerR=$temp_primerR  
  
  #Run cutadapt
  if [[ "${failedMerge_useDirection}" = "NA" ]]; then
    for sample in $(cat ${outdirectory}/sample_order.txt)
    do
        echo "======================================" >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
        echo "On sample: $sample" >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
        echo "======================================" >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
        echo >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
        cutadapt -a "${primerF};required...${revcomp_primerR};optional" \
        -A "${primerR};required...${revcomp_primerF};optional" \
        --discard-untrimmed \
        -j 0 \
        -m 1 \
        -o ${outdirectory}/cutadapt/${sample}_R1_trimmed.fq.gz -p ${outdirectory}/cutadapt/${sample}_R2_trimmed.fq.gz \
        ${readfolderpath}/${sample}_R1.fastq.gz ${readfolderpath}/${sample}_R2.fastq.gz \
        >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt 2>&1
        echo >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
        echo >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
    done
    
    echo "Finished Cutadapt: `date`"
  fi
  
  if [[ "${failedMerge_useDirection}" = "FORWARD" ]]; then
    for sample in $(cat ${outdirectory}/sample_order.txt)
    do
        echo "======================================" >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
        echo "On sample: $sample" >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
        echo "======================================" >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
        echo >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
        cutadapt -a "${primerF};required...${revcomp_primerR};optional" \
        --discard-untrimmed \
        -j 0 \
        -m 1 \
        -o ${outdirectory}/cutadapt/${sample}_R1_trimmed.fq.gz \
        ${readfolderpath}/${sample}_R1.fastq.gz \
        >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt 2>&1
        echo >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
        echo >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
    done
    
    echo "Finished Cutadapt: `date`"
  fi
  
  if [[ "${failedMerge_useDirection}" = "REVERSE" ]]; then
    for sample in $(cat ${outdirectory}/sample_order.txt)
    do
        echo "======================================" >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
        echo "On sample: $sample" >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
        echo "======================================" >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
        echo >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
        cutadapt -a "${primerR};required...${revcomp_primerF};optional" \
        --discard-untrimmed \
        -j 0 \
        -m 1 \
        -o ${outdirectory}/cutadapt/${sample}_R2_trimmed.fq.gz \
        ${readfolderpath}/${sample}_R2.fastq.gz \
        >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt 2>&1
        echo >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
        echo >> ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt
    done
    
    echo "Finished Cutadapt: `date`"
  fi
  
  
  #print cutadapt stats
  echo "Sample	Passing Reads	Passing bp"
  paste ${outdirectory}/sample_order.txt <(grep "passing" ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")") <(grep "filtered" ${outdirectory}/cutadapt/cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")")
  
  if [[ "$bypassflag" = "FALSE" ]]; then
    echo
    echo "Please check Cutadapt success. Proceed? [y/n]"
    read mainmenuinput
    if [[ "$mainmenuinput" = "y" || "$mainmenuinput" = "Y" ]]; then
      echo "Continuing!"
    elif [[ "$mainmenuinput" = "n" || "$mainmenuinput" = "N" ]]; then
      echo "You have chosen to exit"
      exit
    else
      echo "Invalid selection, exiting"
      exit
    fi
  fi

  echo "cutadaptFinished=TRUE" >> ${outdirectory}/progress.txt

fi

#TODO: Deal with samples where ALL reads filtered?

##########################################################################################
##
##    Run DADA2
##
##########################################################################################
if [[ "${dada2_Finished}" = "TRUE" ]]; then
  echo "DADA2 from prior run"
else
  echo
  echo "Running DADA2: `date`"
  if [ -d "${outdirectory}/dada2" ]; then
    rm -r ${outdirectory}/dada2
    mkdir ${outdirectory}/dada2
  else
    mkdir ${outdirectory}/dada2
  fi
  
  echo "Trim and filter in DADA2..."
  
  if [[ "${failedMerge_useDirection}" = "NA" ]]; then
    dada2continue=FALSE
    while [ $dada2continue = FALSE ]
      do
      Rscript --vanilla ${metapipedir}/assets/dada2_step1.R ${workingdirectory}/${outdirectory}/dada2 $dada_minlength $dada_phix $dada_trunQ $dada_maxEE1 $dada_maxEE2 $dada_trimRight $dada_trimLeft \
        1>> ${workingdirectory}/${outdirectory}/dada2/dada2_rscripts_out.log 2>&1
      echo "dada2_step1.R	${workingdirectory}/${outdirectory}/dada2 $dada_minlength $dada_phix $dada_trunQ $dada_maxEE1 $dada_maxEE2 $dada_trimRight $dada_trimLeft" >> ${outdirectory}/Rscript_arguments.log
      echo
      echo "DADA2 Filtering results:"
      echo "Sample	% Reads Passing"
      awk 'NR>1 { print $1, 100 * ( $3 / $2 ) }' ${workingdirectory}/${outdirectory}/dada2/filtered_out_stats.txt
      echo
      echo "Parameters to modify:"
      echo "minLen,rm.phix,truncQ,maxEE-primer1,maxEE-primer2,trimRight,trimLeft"
      echo "Current settings:"
      echo "${dada_minlength},${dada_phix},${dada_trunQ},${dada_maxEE1},${dada_maxEE2},${dada_trimRight},${dada_trimLeft}"
      if [[ "$bypassflag" = "FALSE" ]]; then
        echo "Please check DADA2 filtering success. Proceed? [y/n/m]"
        read mainmenuinput
        if [[ "$mainmenuinput" = "y" || "$mainmenuinput" = "Y" ]]; then
          echo "Continuing!"
          dada2continue=TRUE
        elif [[ "$mainmenuinput" = "m" || "$mainmenuinput" = "M" ]]; then
          echo "You have chosen to redo DADA2 filtering with modified settings."
          echo "Input new settings separated by commas in the same order as above."
          read secondmainmenuinput
          IFS=','
          read -ra ADDR <<< "$secondmainmenuinput"
          dada_minlength=${ADDR[0]}
          dada_phix=${ADDR[1]}
          dada_trunQ=${ADDR[2]}
          dada_maxEE1=${ADDR[3]}
          dada_maxEE2=${ADDR[4]}
          dada_trimRight=${ADDR[5]}
          dada_trimLeft=${ADDR[6]}
          echo "New settings:"
          echo "${dada_minlength},${dada_phix},${dada_trunQ},${dada_maxEE1},${dada_maxEE2},${dada_trimRight},${dada_trimLeft}"
          echo "Rerunning DADA2 filtering"
        elif [[ "$mainmenuinput" = "n" || "$mainmenuinput" = "N" ]]; then
          echo "You have chosen to exit"
          exit
        else
          echo "Invalid selection, try again"
        fi
      else
        dada2continue=TRUE
      fi
    done
  
  
    echo
    echo "Learning error, Dereplication, Merge, and ASVs in DADA2..."
    echo "Please be patient, may take a while. Messages printed to Rscript log."
    echo
    Rscript --vanilla ${metapipedir}/assets/dada2_step2.R ${workingdirectory}/${outdirectory}/dada2 $systemmemoryMB \
        1>> ${workingdirectory}/${outdirectory}/dada2/dada2_rscripts_out.log 2>&1
    echo "dada2_step2.R	${workingdirectory}/${outdirectory}/dada2 $systemmemoryMB" >> ${outdirectory}/Rscript_arguments.log
    
    cat ${outdirectory}/dada2/ASVs_counts.tsv | sed -E 's/^	/x	/' > ${outdirectory}/dada2/ASVs_counts_mod.tsv
    mv ${outdirectory}/dada2/ASVs_counts_mod.tsv ${outdirectory}/dada2/ASVs_counts.tsv
    
    #TODO: minOverlap currently set to 20bp, which seems reasonable. Option for change?
  
    echo
    echo "FINAL DADA2 STATS"
    echo "Note: Please check for a failed merge of forward/reverse sequences"
    echo "Sample	%Reads Retained"
    awk 'NR>1 { print $1, $8 }' ${workingdirectory}/${outdirectory}/dada2/ReadTrimSummary.txt
    if [[ "$bypassflag" = "FALSE" ]]; then
      echo
      echo "Do you wish to Proceed? [y/n]"
      read mainmenuinput
      if [[ "$mainmenuinput" = "y" || "$mainmenuinput" = "Y" ]]; then
        echo "Continuing!"
      elif [[ "$mainmenuinput" = "n" || "$mainmenuinput" = "N" ]]; then
        echo "You have chosen to exit"
        exit
      else
        echo "Invalid selection, exiting"
        exit
      fi
    fi
    
  fi
  
  if [[ "${failedMerge_useDirection}" = "FORWARD" ]]; then
    dada2continue=FALSE
    while [ $dada2continue = FALSE ]
      do
      Rscript --vanilla ${metapipedir}/assets/dada2_step1_mergeFail.R ${workingdirectory}/${outdirectory}/dada2 $dada_minlength $dada_phix $dada_trunQ $dada_maxEE1 $dada_maxEE2 $dada_trimRight $dada_trimLeft forward \
        1>> ${workingdirectory}/${outdirectory}/dada2/dada2_rscripts_out.log 2>&1
      echo "dada2_step1_mergeFail.R	${workingdirectory}/${outdirectory}/dada2 $dada_minlength $dada_phix $dada_trunQ $dada_maxEE1 $dada_maxEE2 $dada_trimRight $dada_trimLeft forward" >> ${outdirectory}/Rscript_arguments.log
      echo
      echo "DADA2 Filtering results:"
      echo "Sample	% Reads Passing"
      awk 'NR>1 { print $1, 100 * ( $3 / $2 ) }' ${workingdirectory}/${outdirectory}/dada2/filtered_out_stats.txt
      echo
      echo "Parameters to modify:"
      echo "minLen,rm.phix,truncQ,maxEE-primer1,maxEE-primer2,trimRight,trimLeft"
      echo "Current settings:"
      echo "${dada_minlength},${dada_phix},${dada_trunQ},${dada_maxEE1},${dada_maxEE2},${dada_trimRight},${dada_trimLeft}"
      if [[ "$bypassflag" = "FALSE" ]]; then
        echo "Please check DADA2 filtering success. Proceed? [y/n/m]"
        read mainmenuinput
        if [[ "$mainmenuinput" = "y" || "$mainmenuinput" = "Y" ]]; then
          echo "Continuing!"
          dada2continue=TRUE
        elif [[ "$mainmenuinput" = "m" || "$mainmenuinput" = "M" ]]; then
          echo "You have chosen to redo DADA2 filtering with modified settings."
          echo "Input new settings separated by commas in the same order as above."
          read secondmainmenuinput
          IFS=','
          read -ra ADDR <<< "$secondmainmenuinput"
          dada_minlength=${ADDR[0]}
          dada_phix=${ADDR[1]}
          dada_trunQ=${ADDR[2]}
          dada_maxEE1=${ADDR[3]}
          dada_maxEE2=${ADDR[4]}
          dada_trimRight=${ADDR[5]}
          dada_trimLeft=${ADDR[6]}
          echo "New settings:"
          echo "${dada_minlength},${dada_phix},${dada_trunQ},${dada_maxEE1},${dada_maxEE2},${dada_trimRight},${dada_trimLeft}"
          echo "Rerunning DADA2 filtering"
        elif [[ "$mainmenuinput" = "n" || "$mainmenuinput" = "N" ]]; then
          echo "You have chosen to exit"
          exit
        else
          echo "Invalid selection, try again"
        fi
      else
        dada2continue=TRUE
      fi
    done
  
  
    echo
    echo "Learning error, Dereplication, and ASVs in DADA2..."
    echo "Please be patient, may take a while. Messages printed to Rscript log."
    echo
    Rscript --vanilla ${metapipedir}/assets/dada2_step2_mergeFail.R ${workingdirectory}/${outdirectory}/dada2 $systemmemoryMB forward \
        1>> ${workingdirectory}/${outdirectory}/dada2/dada2_rscripts_out.log 2>&1
    echo "dada2_step2_mergeFail.R	${workingdirectory}/${outdirectory}/dada2 $systemmemoryMB forward" >> ${outdirectory}/Rscript_arguments.log
    
    cat ${outdirectory}/dada2/ASVs_counts.tsv | sed -E 's/^	/x	/' > ${outdirectory}/dada2/ASVs_counts_mod.tsv
    mv ${outdirectory}/dada2/ASVs_counts_mod.tsv ${outdirectory}/dada2/ASVs_counts.tsv
    
    touch ${outdirectory}/dada2/00_MERGEFAILED_ForwardReadsOnly_usedFor_ASVs_EOM
    
    echo
    echo "FINAL DADA2 STATS"
    echo "Note: Please check for a failed merge of forward/reverse sequences"
    echo "Sample	%Reads Retained"
    awk 'NR>1 { print $1, $6 }' ${workingdirectory}/${outdirectory}/dada2/ReadTrimSummary.txt
    if [[ "$bypassflag" = "FALSE" ]]; then
      echo
      echo "Do you wish to Proceed? [y/n]"
      read mainmenuinput
      if [[ "$mainmenuinput" = "y" || "$mainmenuinput" = "Y" ]]; then
        echo "Continuing!"
      elif [[ "$mainmenuinput" = "n" || "$mainmenuinput" = "N" ]]; then
        echo "You have chosen to exit"
        exit
      else
        echo "Invalid selection, exiting"
        exit
      fi
    fi
    
  fi
  
  if [[ "${failedMerge_useDirection}" = "REVERSE" ]]; then
    dada2continue=FALSE
    while [ $dada2continue = FALSE ]
      do
      Rscript --vanilla ${metapipedir}/assets/dada2_step1_mergeFail.R ${workingdirectory}/${outdirectory}/dada2 $dada_minlength $dada_phix $dada_trunQ $dada_maxEE1 $dada_maxEE2 $dada_trimRight $dada_trimLeft reverse \
        1>> ${workingdirectory}/${outdirectory}/dada2/dada2_rscripts_out.log 2>&1
      echo "dada2_step1_mergeFail.R	${workingdirectory}/${outdirectory}/dada2 $dada_minlength $dada_phix $dada_trunQ $dada_maxEE1 $dada_maxEE2 $dada_trimRight $dada_trimLeft reverse" >> ${outdirectory}/Rscript_arguments.log
      echo
      echo "DADA2 Filtering results:"
      echo "Sample	% Reads Passing"
      awk 'NR>1 { print $1, 100 * ( $3 / $2 ) }' ${workingdirectory}/${outdirectory}/dada2/filtered_out_stats.txt
      echo
      echo "Parameters to modify:"
      echo "minLen,rm.phix,truncQ,maxEE-primer1,maxEE-primer2,trimRight,trimLeft"
      echo "Current settings:"
      echo "${dada_minlength},${dada_phix},${dada_trunQ},${dada_maxEE1},${dada_maxEE2},${dada_trimRight},${dada_trimLeft}"
      if [[ "$bypassflag" = "FALSE" ]]; then
        echo "Please check DADA2 filtering success. Proceed? [y/n/m]"
        read mainmenuinput
        if [[ "$mainmenuinput" = "y" || "$mainmenuinput" = "Y" ]]; then
          echo "Continuing!"
          dada2continue=TRUE
        elif [[ "$mainmenuinput" = "m" || "$mainmenuinput" = "M" ]]; then
          echo "You have chosen to redo DADA2 filtering with modified settings."
          echo "Input new settings separated by commas in the same order as above."
          read secondmainmenuinput
          IFS=','
          read -ra ADDR <<< "$secondmainmenuinput"
          dada_minlength=${ADDR[0]}
          dada_phix=${ADDR[1]}
          dada_trunQ=${ADDR[2]}
          dada_maxEE1=${ADDR[3]}
          dada_maxEE2=${ADDR[4]}
          dada_trimRight=${ADDR[5]}
          dada_trimLeft=${ADDR[6]}
          echo "New settings:"
          echo "${dada_minlength},${dada_phix},${dada_trunQ},${dada_maxEE1},${dada_maxEE2},${dada_trimRight},${dada_trimLeft}"
          echo "Rerunning DADA2 filtering"
        elif [[ "$mainmenuinput" = "n" || "$mainmenuinput" = "N" ]]; then
          echo "You have chosen to exit"
          exit
        else
          echo "Invalid selection, try again"
        fi
      else
        dada2continue=TRUE
      fi
    done
  
  
    echo
    echo "Learning error, Dereplication, and ASVs in DADA2..."
    echo "Please be patient, may take a while. Messages printed to Rscript log."
    echo
    Rscript --vanilla ${metapipedir}/assets/dada2_step2_mergeFail.R ${workingdirectory}/${outdirectory}/dada2 $systemmemoryMB reverse \
        1>> ${workingdirectory}/${outdirectory}/dada2/dada2_rscripts_out.log 2>&1
    echo "dada2_step2_mergeFail.R	${workingdirectory}/${outdirectory}/dada2 $systemmemoryMB reverse" >> ${outdirectory}/Rscript_arguments.log
    
    cat ${outdirectory}/dada2/ASVs_counts.tsv | sed -E 's/^	/x	/' > ${outdirectory}/dada2/ASVs_counts_mod.tsv
    mv ${outdirectory}/dada2/ASVs_counts_mod.tsv ${outdirectory}/dada2/ASVs_counts.tsv
    
    touch ${outdirectory}/dada2/00_MERGEFAILED_ReverseReadsOnly_usedFor_ASVs_EOM
    
    echo
    echo "FINAL DADA2 STATS"
    echo "Note: Please check for a failed merge of forward/reverse sequences"
    echo "Sample	%Reads Retained"
    awk 'NR>1 { print $1, $6 }' ${workingdirectory}/${outdirectory}/dada2/ReadTrimSummary.txt
    if [[ "$bypassflag" = "FALSE" ]]; then
      echo
      echo "Do you wish to Proceed? [y/n]"
      read mainmenuinput
      if [[ "$mainmenuinput" = "y" || "$mainmenuinput" = "Y" ]]; then
        echo "Continuing!"
      elif [[ "$mainmenuinput" = "n" || "$mainmenuinput" = "N" ]]; then
        echo "You have chosen to exit"
        exit
      else
        echo "Invalid selection, exiting"
        exit
      fi
    fi
    
  fi

  echo "dada2_Finished=TRUE" >> ${outdirectory}/progress.txt

fi

##########################################################################################
##
##    Blastn ASVs
##
##########################################################################################
if [[ "${silvaASVflag}" = "TRUE" ]]; then
  echo "Skipping BLASTn to use SILVAngs taxonomy downstream"
else
  if [[ "${blastFinished}" = "TRUE" ]]; then
    echo "BLASTn from prior run"
  else
    if [ -d "${outdirectory}/blast_results" ]; then
      rm -r ${outdirectory}/blast_results
      mkdir ${outdirectory}/blast_results
    else
      mkdir ${outdirectory}/blast_results
    fi
    
    if [[ "$blastflag" = "FALSE" ]]; then
    echo
    echo "Running BLASTn: `date`"
    
    passblastscrutiny=FALSE
    maxtargetseqs=4000
    runthroughcount=3 #temp TODO
    
    while [ $passblastscrutiny = FALSE ]
      do
      let "runthroughcount=runthroughcount+1"
      if [[ "${blastMode}" = "allIN" ]]; then
        blastn -db ${locationNTdatabase}/nt -query ${outdirectory}/dada2/ASVs.fa -outfmt '6 qseqid pident length staxids sacc' -subject_besthit -max_target_seqs $maxtargetseqs -num_threads 6 -out ${outdirectory}/blast_results/ASV_blastn_nt.btab
      elif [[ "${blastMode}" = "mostEnvOUT" ]]; then
        blastn -db ${locationNTdatabase}/nt -query ${outdirectory}/dada2/ASVs.fa -outfmt '6 qseqid pident length staxids sacc' -subject_besthit -max_target_seqs $maxtargetseqs -num_threads 6 -out ${outdirectory}/blast_results/ASV_blastn_nt.btab -negative_taxidlist ${locationNTdatabase}/taxdump/taxid_exclusion_list_leavesinUnclassified.txt
      elif [[ "${blastMode}" = "allEnvOUT" ]]; then
        blastn -db ${locationNTdatabase}/nt -query ${outdirectory}/dada2/ASVs.fa -outfmt '6 qseqid pident length staxids sacc' -subject_besthit -max_target_seqs $maxtargetseqs -num_threads 6 -out ${outdirectory}/blast_results/ASV_blastn_nt.btab -negative_taxidlist ${locationNTdatabase}/taxdump/taxid_exclusion_list_removesUnclassified.txt
      else
        echo "Incorrect blastMode specified"
        exit
      fi
      perl ${metapipedir}/assets/blast_assessment.pl -i ${outdirectory}/blast_results/ASV_blastn_nt.btab -c `grep -c ">" ${outdirectory}/dada2/ASVs.fa` > ${outdirectory}/blast_results/checkmaxtargetseqs.txt
      highestnumber=`cat ${outdirectory}/blast_results/checkmaxtargetseqs.txt | sort -k2 -nr | head -1 | cut -f2`
      if [[ "$highestnumber" -lt "$maxtargetseqs" ]]; then
        passblastscrutiny=TRUE
      elif [[ "$highestnumber" -ge "$maxtargetseqs" && "$runthroughcount" -le 3 ]]; then
        echo "Rerun BLAST Number ${runthroughcount} - Max Target Not High Enough (${maxtargetseqs})"
        let "multiplierseqs=runthroughcount+2"
        let "maxtargetseqs=highestnumber*multiplierseqs"
        echo "New max target seqs = ${maxtargetseqs}"
      else
        echo "Max target seqs inclusive of all top hits cannot be reached"
        echo "Used max target seqs parameter = ${maxtargetseqs}"
        passblastscrutiny=TRUE
      fi
    done
    
    elif [[ "$blastflag" = "TRUE" ]]; then
    echo "Using user supplied BLASTn btab result"
    cp $optionalUserBLASTResult ${outdirectory}/blast_results/ASV_blastn_nt.btab
    fi
  
    echo
    
    echo "blastFinished=TRUE" >> ${outdirectory}/progress.txt
  
  fi
fi

##########################################################################################
##
##    Reformat BLAST output
##
##########################################################################################
if [[ "${silvaASVflag}" = "TRUE" ]]; then
  echo "Skipping BLAST formatting to use SILVAngs taxonomy downstream"
else
  if [[ "${blastformattingFinished}" = "TRUE" ]]; then
    echo "BLAST formatting from prior run"
  else
    if [ -f "${outdirectory}/blast_results/ASV_blastn_nt_formatted.txt" ]; then
      rm ${outdirectory}/blast_results/ASV_blastn_nt_formatted.txt
    fi
  
    echo
    echo "Reformatting BLAST output: `date`"
    
    #let "insertSize=ampliconSize - ${#primerF} - ${#primerR}"
    
    Rscript --vanilla ${metapipedir}/assets/reformat_blast.R ${workingdirectory}/${outdirectory}/blast_results $blastLengthCutoff \
        1>> ${workingdirectory}/${outdirectory}/blast_results/blastreformatting_rscript_out.log 2>&1
    echo "reformat_blast.R	${workingdirectory}/${outdirectory}/blast_results $blastLengthCutoff" >> ${outdirectory}/Rscript_arguments.log
    
    echo
    
    echo "blastformattingFinished=TRUE" >> ${outdirectory}/progress.txt
  
  fi
fi

##########################################################################################
##
##    Running perl script for taxonomy decisions, standard table outs, and KRONA plots
##
##########################################################################################
if [[ "${silvaASVflag}" = "TRUE" ]]; then
  if [[ "${taxonomyscriptFinished}" = "TRUE" ]]; then
    echo "ASV-2-Taxonomy Script results from prior run"
  else
    if [ -d "${outdirectory}/ASV2Taxonomy" ]; then
      rm -r ${outdirectory}/ASV2Taxonomy
      mkdir ${outdirectory}/ASV2Taxonomy
    else
      mkdir ${outdirectory}/ASV2Taxonomy
    fi
    
    echo "Enter the location of the SILVAngs ssu or lsu results directory (i.e. ~/Downloads/results/ssu)"
    read silvaNGSdirectory
    silvangsInputExportFile=${silvaNGSdirectory}/exports/*---otus.csv
    silvangsInputClusterFile=${silvaNGSdirectory}/stats/sequence_cluster_map/data/*.fa.clstr
    echo "Enter the location of the reference taxonomy map for current SILVA database: i.e. tax_slv_ssu_138.1.txt"
    read silvaRefMap
    
    echo
    echo "Running ASV-2-Taxonomy Script: `date`"
    
    cd ${outdirectory}/ASV2Taxonomy
  
    perl ${metapipedir}/assets/asv_taxonomy_processing_figureOuts.pl -a ../dada2/ASVs_counts.tsv \
        -n ${outdirectory} -o ../sample_order.txt -y $silvangsInputExportFile -z $silvangsInputClusterFile -r $silvaRefMap -m $metapipedir -e
    cat ${outdirectory}_asvTaxonomyTable.txt | grep -v "Unknown" > ${outdirectory}_asvTaxonomyTable_NOUNKNOWNS.txt
    cd ${workingdirectory}
    
    cat ${outdirectory}/dada2/ASVs_counts.tsv | grep -v -f ${outdirectory}/ASV2Taxonomy/${outdirectory}_unknown_asvids.txt > ${outdirectory}/ASV2Taxonomy/ASVs_counts_NOUNKNOWNS.tsv
    perl ${metapipedir}/assets/merge_on_taxonomy.pl -a ${outdirectory}/dada2/ASVs_counts.tsv -t ${outdirectory}/ASV2Taxonomy/${outdirectory}_asvTaxonomyTable.txt -o ${outdirectory}/ASV2Taxonomy > ${outdirectory}/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy.tsv
    cat ${outdirectory}/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy.tsv | grep -v -f ${outdirectory}/ASV2Taxonomy/${outdirectory}_unknown_asvids.txt > ${outdirectory}/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy_NOUNKNOWNS.tsv
    mkdir ${outdirectory}/ASV2Taxonomy/KRONA_plots
    mkdir ${outdirectory}/ASV2Taxonomy/KRONA_plots/KRONA_inputs
    mv ${outdirectory}/ASV2Taxonomy/MP* ${outdirectory}/ASV2Taxonomy/KRONA_plots/KRONA_inputs/
    mv ${outdirectory}/ASV2Taxonomy/${outdirectory}_wholeKRONA.txt ${outdirectory}/ASV2Taxonomy/KRONA_plots/KRONA_inputs/${outdirectory}_samplesSummedKRONA.txt
    mv ${outdirectory}/ASV2Taxonomy/${outdirectory}_master_krona.html ${outdirectory}/ASV2Taxonomy/KRONA_plots/
    mv ${outdirectory}/ASV2Taxonomy/${outdirectory}_wholeKRONA.html ${outdirectory}/ASV2Taxonomy/KRONA_plots/${outdirectory}_samplesSummedKRONA.html
    mv ${outdirectory}/ASV2Taxonomy/${outdirectory}_allin_KRONA.txt ${outdirectory}/ASV2Taxonomy/${outdirectory}_allin_TaxonomyASVSampleCount_byline.txt
    
    perl ${metapipedir}/assets/stats.pl -a ${outdirectory}/dada2/ASVs_counts.tsv -t ${outdirectory}/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy.tsv -i ${outdirectory}/ASV2Taxonomy/${outdirectory}_asvTaxonomyTable.txt > ${outdirectory}/ASV2Taxonomy/basic_ASV_taxonomy_stats.txt
    
    echo "taxonomyscriptFinished=TRUE" >> ${outdirectory}/progress.txt
  
  fi
  
else
  if [[ "${taxonomyscriptFinished}" = "TRUE" ]]; then
    echo "ASV-2-Taxonomy Script results from prior run"
  else
    if [ -d "${outdirectory}/ASV2Taxonomy" ]; then
      rm -r ${outdirectory}/ASV2Taxonomy
      mkdir ${outdirectory}/ASV2Taxonomy
      if [[ "$removeASVsFILE" != "NULL" ]]; then
        mv ${outdirectory}/dada2/ASVs_unfiltered.fa ${outdirectory}/dada2/ASVs.fa
        mv ${outdirectory}/dada2/ASVs_counts_unfiltered.tsv ${outdirectory}/dada2/ASVs_counts.tsv
        #WILL THROW AN ERROR: If you run with removeASVsFILE given then remove that parameter without moving the *unfiltered* back to their original names.
      fi
    else
      mkdir ${outdirectory}/ASV2Taxonomy
    fi
  
    echo
    echo "Running ASV-2-Taxonomy Script: `date`"
    
    #prep taxon ids for taxonkit
    cat ${outdirectory}/blast_results/ASV_blastn_nt_formatted.txt | cut -f4 | sed -E "s/ //g" | tr ',' '\n' | tr ';' '\n' | sort | uniq | grep -v "taxid" > ${outdirectory}/ASV2Taxonomy/taxids.txt
    taxonkit lineage ${outdirectory}/ASV2Taxonomy/taxids.txt | awk '$2!=""' > ${outdirectory}/ASV2Taxonomy/taxonkit_out.txt
    taxonkit reformat ${outdirectory}/ASV2Taxonomy/taxonkit_out.txt | cut -f1,3 > ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt
    
    if [[ "$bypassflag" = "FALSE" ]]; then
      echo
      echo "Reformatted taxon strings created. Options:"
      echo "Continue without changes [c]"
      echo "Manually edit file and replace in same location with identical file structure [m]"
      echo "    (Make choice when file is modified and you are ready to proceed)"
      echo "Automatically fill gaps in reformatted taxonkit hierarchy [a]"
      read mainmenuinput
      if [[ "$mainmenuinput" = "c" || "$mainmenuinput" = "C" ]]; then
        cat ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt | sed -E 's/[^A-Za-z0-9;[:blank:]]/_/g' > ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt_temp
        mv ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out_ORIGINAL.txt
        mv ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt_temp ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt
        echo "Continuing!"
      elif [[ "$mainmenuinput" = "m" || "$mainmenuinput" = "M" ]]; then
        echo "Continuing!"
        cat ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt | tr '\r' '\n' | tr -s '\n' > ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt_temp
        cat ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt_temp | sed -E 's/[^A-Za-z0-9;[:blank:]]/_/g' > ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt
        rm ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt_temp
      elif [[ "$mainmenuinput" = "a" || "$mainmenuinput" = "A" ]]; then
        echo "Reformatting..."
        echo "Original reformatted taxonkit out stored at ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out_ORIGINAL.txt"
        perl ${metapipedir}/assets/fillIn_taxonkit.pl -i ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt > ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt_temp
        mv ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out_ORIGINAL.txt
        cat ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt_temp | sed -E 's/[^A-Za-z0-9;[:blank:]]/_/g' > ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt
        rm ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt_temp
        echo "Continuing!"
      else
        echo "Invalid selection, exiting"
        exit
      fi
    elif [[ "$bypassflag" = "TRUE" ]]; then
      echo "Automatically reformatting taxonkit hierarchy..."
      echo "Original reformatted taxonkit out stored at ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out_ORIGINAL.txt"
      perl ${metapipedir}/assets/fillIn_taxonkit.pl -i ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt > ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt_temp
      mv ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out_ORIGINAL.txt
      cat ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt_temp | sed -E 's/[^A-Za-z0-9;[:blank:]]/_/g' > ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt
      rm ${outdirectory}/ASV2Taxonomy/reformatted_taxonkit_out.txt_temp
      echo "Continuing!"
    fi
      
    
    if [[ "$removeASVsFILE" = "NULL" ]]; then
      cd ${outdirectory}/ASV2Taxonomy
    
      perl ${metapipedir}/assets/asv_taxonomy_processing_figureOuts.pl -a ../dada2/ASVs_counts.tsv -s ../blast_results/ASV_blastn_nt_formatted.txt \
          -t reformatted_taxonkit_out.txt -f $taxonomyCutoffs -n ${outdirectory} -c ${locationNTdatabase}/taxdump/common_names.dmp -o ../sample_order.txt
      cat ${outdirectory}_asvTaxonomyTable.txt | grep -v "Unknown" > ${outdirectory}_asvTaxonomyTable_NOUNKNOWNS.txt
      cd ${workingdirectory}
      cat ${outdirectory}/ASV2Taxonomy/${outdirectory}_unknown_asvids.txt | cut -f1 | sed -E 's/$/	/' > ${outdirectory}/ASV2Taxonomy/temp_grep_unknowns
      cat ${outdirectory}/dada2/ASVs_counts.tsv | grep -v -f ${outdirectory}/ASV2Taxonomy/temp_grep_unknowns > ${outdirectory}/ASV2Taxonomy/ASVs_counts_NOUNKNOWNS.tsv
      perl ${metapipedir}/assets/merge_on_taxonomy.pl -a ${outdirectory}/dada2/ASVs_counts.tsv -t ${outdirectory}/ASV2Taxonomy/${outdirectory}_asvTaxonomyTable.txt -o ${outdirectory}/ASV2Taxonomy > ${outdirectory}/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy.tsv
      cat ${outdirectory}/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy.tsv | grep -v -f ${outdirectory}/ASV2Taxonomy/temp_grep_unknowns > ${outdirectory}/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy_NOUNKNOWNS.tsv
      rm ${outdirectory}/ASV2Taxonomy/temp_grep_unknowns
      mkdir ${outdirectory}/ASV2Taxonomy/KRONA_plots
      mkdir ${outdirectory}/ASV2Taxonomy/KRONA_plots/KRONA_inputs
      mv ${outdirectory}/ASV2Taxonomy/MP* ${outdirectory}/ASV2Taxonomy/KRONA_plots/KRONA_inputs/
      mv ${outdirectory}/ASV2Taxonomy/${outdirectory}_wholeKRONA.txt ${outdirectory}/ASV2Taxonomy/KRONA_plots/KRONA_inputs/${outdirectory}_samplesSummedKRONA.txt
      mv ${outdirectory}/ASV2Taxonomy/${outdirectory}_master_krona.html ${outdirectory}/ASV2Taxonomy/KRONA_plots/
      mv ${outdirectory}/ASV2Taxonomy/${outdirectory}_wholeKRONA.html ${outdirectory}/ASV2Taxonomy/KRONA_plots/${outdirectory}_samplesSummedKRONA.html
      mv ${outdirectory}/ASV2Taxonomy/${outdirectory}_allin_KRONA.txt ${outdirectory}/ASV2Taxonomy/${outdirectory}_allin_TaxonomyASVSampleCount_byline.txt
      
      perl ${metapipedir}/assets/stats.pl -a ${outdirectory}/dada2/ASVs_counts.tsv -t ${outdirectory}/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy.tsv -i ${outdirectory}/ASV2Taxonomy/${outdirectory}_asvTaxonomyTable.txt > ${outdirectory}/ASV2Taxonomy/basic_ASV_taxonomy_stats.txt
      
    else
      cd ${outdirectory}/ASV2Taxonomy
      
      perl ${metapipedir}/assets/asv_taxonomy_processing_figureOuts.pl -a ../dada2/ASVs_counts.tsv -s ../blast_results/ASV_blastn_nt_formatted.txt -d $removeASVsFILE \
          -t reformatted_taxonkit_out.txt -f $taxonomyCutoffs -n ${outdirectory} -c ${locationNTdatabase}/taxdump/common_names.dmp -o ../sample_order.txt
      
      #Move non-filtered files to folder, replace with files from "IGNORING_ASVs" directory
      mkdir unfiltered_files
      mkdir -p KRONA_plots/KRONA_inputs
      mv MP* KRONA_plots/KRONA_inputs/
      mv ${outdirectory}_wholeKRONA.txt KRONA_plots/KRONA_inputs/${outdirectory}_samplesSummedKRONA.txt
      mv ${outdirectory}_master_krona.html KRONA_plots/
      mv ${outdirectory}_wholeKRONA.html KRONA_plots/${outdirectory}_samplesSummedKRONA.html
      mv ${outdirectory}_allin_KRONA.txt ${outdirectory}_allin_TaxonomyASVSampleCount_byline.txt
      mv KRONA_plots unfiltered_files/
      mv ${outdirectory}_allin_TaxonomyASVSampleCount_byline.txt unfiltered_files/
      mv ${outdirectory}_asvTaxonomyTable.txt unfiltered_files/
      mv ${outdirectory}_barchart_forR.txt unfiltered_files/
      mv ${outdirectory}_barchart.txt unfiltered_files/
      mv ${outdirectory}_heatmap_multiASV.txt unfiltered_files/
      mv ${outdirectory}_NO_UNKNOWNS_barchart.txt unfiltered_files/
      
      mkdir -p IGNORING_ASVs/KRONA_plots/KRONA_inputs
      mv IGNORING_ASVs/MP* IGNORING_ASVs/KRONA_plots/KRONA_inputs/
      mv IGNORING_ASVs/${outdirectory}_IGNORE_master_krona.html IGNORING_ASVs/KRONA_plots/
      mv IGNORING_ASVs/${outdirectory}_IGNORE_allin_KRONA.txt ${outdirectory}_allin_TaxonomyASVSampleCount_byline.txt
      mv IGNORING_ASVs/KRONA_plots ./
      mv IGNORING_ASVs/${outdirectory}_ASVs_to_IGNORE.txt ./
      mv IGNORING_ASVs/${outdirectory}_IGNORE_asvTaxonomyTable.txt ./${outdirectory}_asvTaxonomyTable.txt
      mv IGNORING_ASVs/${outdirectory}_IGNORE_barchart_forR.txt ./${outdirectory}_barchart_forR.txt
      mv IGNORING_ASVs/${outdirectory}_IGNORE_barchart.txt ./${outdirectory}_barchart.txt
      mv IGNORING_ASVs/${outdirectory}_IGNORE_heatmap_multiASV.txt ./${outdirectory}_heatmap_multiASV.txt
      mv IGNORING_ASVs/${outdirectory}_IGNORE_NO_UNKNOWNS_barchart.txt ./${outdirectory}_NO_UNKNOWNS_barchart.txt
      rmdir IGNORING_ASVs
      
      cat ${outdirectory}_asvTaxonomyTable.txt | grep -v "Unknown" > ${outdirectory}_asvTaxonomyTable_NOUNKNOWNS.txt
      cd ${workingdirectory}
      
      mv ${outdirectory}/dada2/ASVs.fa ${outdirectory}/dada2/ASVs_unfiltered.fa
      mv ${outdirectory}/dada2/ASVs_counts.tsv ${outdirectory}/dada2/ASVs_counts_unfiltered.tsv
      
      cat ${outdirectory}/ASV2Taxonomy/${outdirectory}_ASVs_to_IGNORE.txt | sed -E 's/$/\$/' > ${outdirectory}/temp
      grep -A 1 -f ${outdirectory}/temp ${outdirectory}/dada2/ASVs_unfiltered.fa | sed -E 's/$/\$/' > ${outdirectory}/temp2
      grep -v "^--\$$" ${outdirectory}/temp2 > ${outdirectory}/temp
      grep -v -f ${outdirectory}/temp ${outdirectory}/dada2/ASVs_unfiltered.fa > ${outdirectory}/dada2/ASVs.fa
      rm ${outdirectory}/temp
      rm ${outdirectory}/temp2
      #Note for potential bug: Will remove the sequence (but not the header) for sequences identical to an ASV in the IGNORE list. 
      
      cat ${outdirectory}/ASV2Taxonomy/${outdirectory}_ASVs_to_IGNORE.txt | sed -E 's/$/	/' > ${outdirectory}/temp
      cat ${outdirectory}/dada2/ASVs_counts_unfiltered.tsv | grep -v -f ${outdirectory}/temp > ${outdirectory}/dada2/ASVs_counts.tsv
      rm ${outdirectory}/temp
      
      cat ${outdirectory}/ASV2Taxonomy/${outdirectory}_unknown_asvids.txt | cut -f1 | sed -E 's/$/	/' > ${outdirectory}/ASV2Taxonomy/temp_grep_unknowns
      cat ${outdirectory}/dada2/ASVs_counts.tsv | grep -v -f ${outdirectory}/ASV2Taxonomy/temp_grep_unknowns > ${outdirectory}/ASV2Taxonomy/ASVs_counts_NOUNKNOWNS.tsv
      perl ${metapipedir}/assets/merge_on_taxonomy.pl -a ${outdirectory}/dada2/ASVs_counts.tsv -t ${outdirectory}/ASV2Taxonomy/${outdirectory}_asvTaxonomyTable.txt -o ${outdirectory}/ASV2Taxonomy > ${outdirectory}/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy.tsv
      cat ${outdirectory}/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy.tsv | grep -v -f ${outdirectory}/ASV2Taxonomy/temp_grep_unknowns > ${outdirectory}/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy_NOUNKNOWNS.tsv
      rm ${outdirectory}/ASV2Taxonomy/temp_grep_unknowns
      
      perl ${metapipedir}/assets/stats.pl -a ${outdirectory}/dada2/ASVs_counts.tsv -t ${outdirectory}/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy.tsv -i ${outdirectory}/ASV2Taxonomy/${outdirectory}_asvTaxonomyTable.txt > ${outdirectory}/ASV2Taxonomy/basic_ASV_taxonomy_stats.txt
      
    fi
  
    echo "taxonomyscriptFinished=TRUE" >> ${outdirectory}/progress.txt
  
  fi
fi
##########################################################################################
##
##    File cleanup
##
##########################################################################################
if [[ "$keepIntermediateFiles" = "FALSE" ]]; then
  rm -f ${outdirectory}/cutadapt/*_trimmed.fq.gz
  rm -f ${outdirectory}/dada2/*_filtered.fq.gz
  if [[ "${silvaASVflag}" != "TRUE" ]]; then
    gzip -q -9 ${outdirectory}/blast_results/ASV_blastn_nt.btab
  fi
fi

##########################################################################################
##
##    FIGURES & ANALYSIS
##
##########################################################################################
#Defaults, user can override:
filterPercent=5
filterLowQualSamples="FALSE"
filterPercentLowQualSamples=30
pieScale=6 #TODO assign appropriate default
removeNA="FALSE"
providedTaxaOfInterest="FALSE"
taxaOfInterestFile="NULL"
taxaOfInterestLevel="NULL"
#Interpreted from sample metadata file:
controlPos="FALSE"
controlNeg="FALSE"
replicates="FALSE"
sites="FALSE"
chemData="FALSE"
controlspresent="FALSE"
groupsDefinedFlag="FALSE"
locationChemHeaders="NULL"

source $figureparamfilepath

cp $figureparamfilepath ${outdirectory}/figure_config_file.txt

#TODO Check previous and current config files and skip this if finished. Likely will want this to be more granular.

if [[ "${figuresFinished}" = "TRUE" ]]; then
  echo "Figure creation from prior run"
else
  if [ -d "${outdirectory}/Figures" ]; then
    rm -r ${outdirectory}/Figures
    mkdir ${outdirectory}/Figures
  else
    mkdir ${outdirectory}/Figures
  fi
  if [ -d "${outdirectory}/processed_tables" ]; then
    rm -r ${outdirectory}/processed_tables
    mkdir ${outdirectory}/processed_tables
  else
    mkdir ${outdirectory}/processed_tables
  fi

cp -r ${outdirectory}/ASV2Taxonomy/KRONA_plots ${outdirectory}/Figures/00_KRONA_plots
rm -r ${outdirectory}/Figures/00_KRONA_plots/KRONA_inputs
  
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
      elif [[ "$header" = "controls" ]]; then
      positiveCount=`cat ${workingdirectory}/${outdirectory}/temp | grep -c "positive"`
      negativeCount=`cat ${workingdirectory}/${outdirectory}/temp | grep -c "negative"`
      if [[ $positiveCount -ge 1 ]]; then
        controlPos="TRUE"
      fi
      if [[ $negativeCount -ge 1 ]]; then
        controlNeg="TRUE"
      fi
      if [[ $positiveCount -ge 1 || $negativeCount -ge 1 ]]; then
        controlspresent="TRUE"
      fi
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

#echo "Replicates: $replicates"
#echo "Sites: $sites"
#echo "Controls: $controlspresent"
#echo "Positive Control: $controlPos"
#echo "Negative Control: $controlNeg"
#echo "Groups called: $groupsDefinedFlag"
#echo "Number of groups: $numberGroupsDefined"
#echo "Chem data: $chemData"
#echo "Location: $locationChemHeaders"

#Maps
mkdir ${outdirectory}/Figures/01_Maps
Rscript --vanilla ${metapipedir}/assets/maps.R ${workingdirectory}/${outdirectory}/Figures/01_Maps ${workingdirectory}/${outdirectory}/sample_metadata_forR.txt $replicates $sites ${workingdirectory}/${outdirectory}/ASV2Taxonomy/${outdirectory}_NO_UNKNOWNS_barchart.txt $filterPercent $pieScale \
    1>> ${workingdirectory}/${outdirectory}/Figures/01_Maps/maps_rscript_out.log 2>&1
echo "maps.R	${workingdirectory}/${outdirectory}/Figures/01_Maps ${workingdirectory}/${outdirectory}/sample_metadata_forR.txt $replicates $sites ${workingdirectory}/${outdirectory}/ASV2Taxonomy/${outdirectory}_NO_UNKNOWNS_barchart.txt $filterPercent $pieScale" >> ${outdirectory}/Rscript_arguments.log

rm -f ${workingdirectory}/${outdirectory}/Figures/01_Maps/Rplot*

#Tables
perl ${metapipedir}/assets/barchart_filterLowAbund.pl -i ${outdirectory}/ASV2Taxonomy/${outdirectory}_barchart_forR.txt -f $filterPercent > ${outdirectory}/ASV2Taxonomy/${outdirectory}_barchart_forR_filtLowAbund_zzOther.txt

Rscript --vanilla ${metapipedir}/assets/process_tables.R ${workingdirectory}/${outdirectory} ${workingdirectory}/${outdirectory}/ASV2Taxonomy/ASVs_counts_NOUNKNOWNS.tsv ${workingdirectory}/${outdirectory}/ASV2Taxonomy/${outdirectory}_asvTaxonomyTable_NOUNKNOWNS.txt ${workingdirectory}/${outdirectory}/sample_metadata_forR.txt $filterPercent $controlspresent $filterLowQualSamples $filterPercentLowQualSamples ${workingdirectory}/${outdirectory}/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy_NOUNKNOWNS.tsv $sites $replicates ${workingdirectory}/${outdirectory}/dada2/ASVs_counts.tsv \
  1>> ${workingdirectory}/${outdirectory}/processed_tables/table_rscript_out.log 2>&1
echo "process_tables.R	${workingdirectory}/${outdirectory} ${workingdirectory}/${outdirectory}/ASV2Taxonomy/ASVs_counts_NOUNKNOWNS.tsv ${workingdirectory}/${outdirectory}/ASV2Taxonomy/${outdirectory}_asvTaxonomyTable_NOUNKNOWNS.txt ${workingdirectory}/${outdirectory}/sample_metadata_forR.txt $filterPercent $controlspresent $filterLowQualSamples $filterPercentLowQualSamples ${workingdirectory}/${outdirectory}/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy_NOUNKNOWNS.tsv $sites $replicates ${workingdirectory}/${outdirectory}/dada2/ASVs_counts.tsv" >> ${outdirectory}/Rscript_arguments.log
perl ${metapipedir}/assets/filter_lowabundance_taxa.pl -a ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.tsv -t ${outdirectory}/ASV2Taxonomy/${outdirectory}_asvTaxonomyTable_NOUNKNOWNS.txt -p $filterPercent > ${outdirectory}/processed_tables/ASVTaxonomyTable_NOUNKNOWNS_replaceLowAbund2zzOther.txt

if [[ "${controlspresent}" = "TRUE" ]]; then
    if [[ "${filterLowQualSamples}" = "TRUE" ]]; then
      compareproc1=`cat ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_controlsRemoved.tsv`
      compareproc2=`cat ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_lowEffortSamplesRemoved.tsv`
      compareproc3=`cat ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_controlsRemoved.tsv`
      compareproc4=`cat ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_lowEffortSamplesRemoved.tsv`
      compareproc5=`cat ${outdirectory}/processed_tables/ASVs_counts_lowEffortSamplesRemoved.tsv`
      compareproc6=`cat ${outdirectory}/processed_tables/ASVs_counts_controlsRemoved.tsv`
      
    
      if [[ "${compareproc1}" = "${compareproc2}" ]]; then
        echo "Low quality sample filtering ON, but none met criteria to be removed"
        filterLowQualSamples=FALSE
        rm ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_lowEffortSamplesRemoved.tsv
      fi
      if [[ "${compareproc3}" = "${compareproc4}" ]]; then
        rm ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_lowEffortSamplesRemoved.tsv
      fi
      if [[ "${compareproc5}" = "${compareproc6}" ]]; then
        rm ${outdirectory}/processed_tables/ASVs_counts_lowEffortSamplesRemoved.tsv
      fi
    fi
else #CONTROLSPRESENT = FALSE
    if [[ "${filterLowQualSamples}" = "TRUE" ]]; then
      compareproc1=`cat ${outdirectory}/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy_NOUNKNOWNS.tsv`
      compareproc2=`cat ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_lowEffortSamplesRemoved.tsv`
      compareproc3=`cat ${outdirectory}/ASV2Taxonomy/ASVs_counts_NOUNKNOWNS.tsv`
      compareproc4=`cat ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_lowEffortSamplesRemoved.tsv`
      compareproc5=`cat ${outdirectory}/dada2/ASVs_counts.tsv`
      compareproc6=`cat ${outdirectory}/processed_tables/ASVs_counts_lowEffortSamplesRemoved.tsv`
      
      if [[ "${compareproc1}" = "${compareproc2}" ]]; then
        echo "Low quality sample filtering ON, but none met criteria to be removed"
        filterLowQualSamples=FALSE
        rm ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_lowEffortSamplesRemoved.tsv
      fi
      if [[ "${compareproc3}" = "${compareproc4}" ]]; then
        rm ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_lowEffortSamplesRemoved.tsv
      fi
      if [[ "${compareproc5}" = "${compareproc6}" ]]; then
        rm ${outdirectory}/processed_tables/ASVs_counts_lowEffortSamplesRemoved.tsv
      fi
    fi
fi

for f in ${workingdirectory}/${outdirectory}/processed_tables/sample_metadata*; do cat $f | sed -E 's/, /_/g' > ${f}_mod; mv ${f}_mod $f; done

#Phyloseq figures
mkdir -p ${outdirectory}/Figures/02_Barcharts/read_count
mkdir -p ${outdirectory}/Figures/02_Barcharts/relative_abundance
mkdir -p ${outdirectory}/Figures/03_Heatmaps/ASV_based
mkdir -p ${outdirectory}/Figures/03_Heatmaps/Taxonomy_merge_based
mkdir -p ${outdirectory}/Figures/04_Alpha_Diversity/ASV_based
mkdir -p ${outdirectory}/Figures/04_Alpha_Diversity/Taxonomy_merge_based
mkdir -p ${outdirectory}/Figures/05_Ordination/ASV_based/read_count
mkdir -p ${outdirectory}/Figures/05_Ordination/ASV_based/relative_abundance
mkdir -p ${outdirectory}/Figures/05_Ordination/Taxonomy_merge_based/read_count
mkdir -p ${outdirectory}/Figures/05_Ordination/Taxonomy_merge_based/relative_abundance
mkdir -p ${outdirectory}/Figures/05_Ordination/Taxonomy_merge_based/filterInclude_TOSPECIES_only/read_count
mkdir -p ${outdirectory}/Figures/05_Ordination/Taxonomy_merge_based/filterInclude_TOSPECIES_only/relative_abundance
mkdir -p ${outdirectory}/Figures/06_Network/ASV_based/read_count
mkdir -p ${outdirectory}/Figures/06_Network/ASV_based/relative_abundance
mkdir -p ${outdirectory}/Figures/06_Network/Taxonomy_merge_based/read_count
mkdir -p ${outdirectory}/Figures/06_Network/Taxonomy_merge_based/relative_abundance
mkdir -p ${outdirectory}/Figures/07_Rarefaction_Curves
if [[ "${providedTaxaOfInterest}" = "TRUE" ]]; then
mkdir -p ${outdirectory}/Figures/Taxa_of_interest/02_Barcharts/read_count
mkdir -p ${outdirectory}/Figures/Taxa_of_interest/02_Barcharts/relative_abundance
mkdir -p ${outdirectory}/Figures/Taxa_of_interest/03_Heatmaps/ASV_based
mkdir -p ${outdirectory}/Figures/Taxa_of_interest/03_Heatmaps/Taxonomy_merge_based
mkdir -p ${outdirectory}/Figures/Taxa_of_interest/06_Network/ASV_based/read_count
mkdir -p ${outdirectory}/Figures/Taxa_of_interest/06_Network/ASV_based/relative_abundance
mkdir -p ${outdirectory}/Figures/Taxa_of_interest/06_Network/Taxonomy_merge_based/read_count
mkdir -p ${outdirectory}/Figures/Taxa_of_interest/06_Network/Taxonomy_merge_based/relative_abundance
fi

#echo ${workingdirectory}/${outdirectory}/Figures ${workingdirectory} ${outdirectory} $controlspresent $filterLowQualSamples $replicates $sites $filterPercent $removeNA $providedTaxaOfInterest $groupsDefinedFlag $numberGroupsDefined $taxaOfInterestLevel $taxaOfInterestFile $chemData $locationChemHeaders
Rscript --vanilla ${metapipedir}/assets/phyloseq.R ${workingdirectory}/${outdirectory}/Figures ${workingdirectory} ${outdirectory} $controlspresent $filterLowQualSamples $replicates $sites $filterPercent $removeNA $providedTaxaOfInterest $groupsDefinedFlag $numberGroupsDefined $taxaOfInterestLevel $taxaOfInterestFile $chemData $locationChemHeaders \
  1>> ${workingdirectory}/${outdirectory}/Figures/phyloseq_rscript_out.log 2>&1
echo "phyloseq.R	${workingdirectory}/${outdirectory}/Figures ${workingdirectory} ${outdirectory} $controlspresent $filterLowQualSamples $replicates $sites $filterPercent $removeNA $providedTaxaOfInterest $groupsDefinedFlag $numberGroupsDefined $taxaOfInterestLevel $taxaOfInterestFile $chemData $locationChemHeaders" >> ${outdirectory}/Rscript_arguments.log

rm -f ${workingdirectory}/${outdirectory}/Figures/02_Barcharts/read_count/Rplots.pdf

Rscript --vanilla ${metapipedir}/assets/barchart_terminaltaxa.R ${workingdirectory}/${outdirectory}/Figures ${workingdirectory}/${outdirectory}/ASV2Taxonomy/${outdirectory}_barchart_forR.txt ${workingdirectory}/${outdirectory}/sample_order.txt ${workingdirectory}/${outdirectory}/ASV2Taxonomy/${outdirectory}_barchart_forR_filtLowAbund_zzOther.txt \
 1>> ${workingdirectory}/${outdirectory}/Figures/phyloseq_rscript_out.log 2>&1
echo "barchart_terminaltaxa.R	${workingdirectory}/${outdirectory}/Figures ${workingdirectory}/${outdirectory}/ASV2Taxonomy/${outdirectory}_barchart_forR.txt ${workingdirectory}/${outdirectory}/sample_order.txt ${workingdirectory}/${outdirectory}/ASV2Taxonomy/${outdirectory}_barchart_forR_filtLowAbund_zzOther.txt" >> ${outdirectory}/Rscript_arguments.log

rm -f ${workingdirectory}/${outdirectory}/Figures/02_Barcharts/relative_abundance/Rplots.pdf

mkdir -p ${outdirectory}/Figures/08_EnvironmentFit_Ordination/ASV_based
mkdir -p ${outdirectory}/Figures/08_EnvironmentFit_Ordination/Taxonomy_merge_based
Rscript --vanilla ${metapipedir}/assets/environment_fit_ordination.R ${workingdirectory}/${outdirectory}/Figures ${workingdirectory}/${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_percentabund.tsv ${workingdirectory}/${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.tsv ${workingdirectory}/${outdirectory}/sample_metadata_forR.txt $replicates $sites $chemData $locationChemHeaders \
  1>> ${workingdirectory}/${outdirectory}/Figures/08_EnvironmentFit_Ordination/envfit_rscript_out.log 2>&1
echo "environment_fit_ordination.R	${workingdirectory}/${outdirectory}/Figures ${workingdirectory}/${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_percentabund.tsv ${workingdirectory}/${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_percentabund.tsv ${workingdirectory}/${outdirectory}/sample_metadata_forR.txt $replicates $sites $chemData $locationChemHeaders" >> ${outdirectory}/Rscript_arguments.log

#Replicate section
if [[ "${replicates}" = "TRUE" ]]; then
  mkdir ${outdirectory}/processed_tables/replicate_based_detection
  if [[ "${controlspresent}" = "TRUE" ]]; then
    if [[ "${filterLowQualSamples}" = "TRUE" ]]; then
      perl ${metapipedir}/assets/replicate_presence_absence.pl -i ${outdirectory}/processed_tables/ASVs_counts_lowEffortSamplesRemoved.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/presenceabsence_unmarked_ASVbased_withUnknowns_filtsamples.txt
      perl ${metapipedir}/assets/replicate_presence_absence.pl -i ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_lowEffortSamplesRemoved.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/presenceabsence_unmarked_ASVbased_NoUnknowns_filtsamples.txt
      perl ${metapipedir}/assets/replicate_presence_absence.pl -i ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_lowEffortSamplesRemoved.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/presenceabsence_unmarked_TAXAbased_NoUnknowns_filtsamples.txt
      perl ${metapipedir}/assets/replicate_compRelabund_detection.pl -i ${outdirectory}/processed_tables/ASVs_counts_lowEffortSamplesRemoved.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_ASVbased_withUnknowns_filtsamples.txt
      perl ${metapipedir}/assets/replicate_compRelabund_detection.pl -i ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_lowEffortSamplesRemoved.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_ASVbased_NoUnknowns_filtsamples.txt
      perl ${metapipedir}/assets/replicate_compRelabund_detection.pl -i ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_lowEffortSamplesRemoved.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_TAXAbased_NoUnknowns_filtsamples.txt

    fi
    if [[ "${filterLowQualSamples}" = "FALSE" ]]; then
      perl ${metapipedir}/assets/replicate_presence_absence.pl -i ${outdirectory}/processed_tables/ASVs_counts_controlsRemoved.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/presenceabsence_unmarked_ASVbased_withUnknowns_filtsamples.txt
      perl ${metapipedir}/assets/replicate_presence_absence.pl -i ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_controlsRemoved.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/presenceabsence_unmarked_ASVbased_NoUnknowns_filtsamples.txt
      perl ${metapipedir}/assets/replicate_presence_absence.pl -i ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_controlsRemoved.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/presenceabsence_unmarked_TAXAbased_NoUnknowns_filtsamples.txt
      perl ${metapipedir}/assets/replicate_compRelabund_detection.pl -i ${outdirectory}/processed_tables/ASVs_counts_controlsRemoved.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_ASVbased_withUnknowns_filtsamples.txt
      perl ${metapipedir}/assets/replicate_compRelabund_detection.pl -i ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_controlsRemoved.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_ASVbased_NoUnknowns_filtsamples.txt
      perl ${metapipedir}/assets/replicate_compRelabund_detection.pl -i ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_controlsRemoved.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_TAXAbased_NoUnknowns_filtsamples.txt
    fi
  fi
  if [[ "${controlspresent}" = "FALSE" ]]; then
    if [[ "${filterLowQualSamples}" = "TRUE" ]]; then
      perl ${metapipedir}/assets/replicate_presence_absence.pl -i ${outdirectory}/processed_tables/ASVs_counts_lowEffortSamplesRemoved.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/presenceabsence_unmarked_ASVbased_withUnknowns_filtsamples.txt
      perl ${metapipedir}/assets/replicate_presence_absence.pl -i ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_lowEffortSamplesRemoved.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/presenceabsence_unmarked_ASVbased_NoUnknowns_filtsamples.txt
      perl ${metapipedir}/assets/replicate_presence_absence.pl -i ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_lowEffortSamplesRemoved.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/presenceabsence_unmarked_TAXAbased_NoUnknowns_filtsamples.txt
      perl ${metapipedir}/assets/replicate_compRelabund_detection.pl -i ${outdirectory}/processed_tables/ASVs_counts_lowEffortSamplesRemoved.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_ASVbased_withUnknowns_filtsamples.txt
      perl ${metapipedir}/assets/replicate_compRelabund_detection.pl -i ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_lowEffortSamplesRemoved.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_ASVbased_NoUnknowns_filtsamples.txt
      perl ${metapipedir}/assets/replicate_compRelabund_detection.pl -i ${outdirectory}/processed_tables/ASVs_counts_NOUNKNOWNS_collapsedOnTaxonomy_lowEffortSamplesRemoved.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_TAXAbased_NoUnknowns_filtsamples.txt
    fi
    if [[ "${filterLowQualSamples}" = "FALSE" ]]; then
      perl ${metapipedir}/assets/replicate_presence_absence.pl -i ${outdirectory}/dada2/ASVs_counts.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/presenceabsence_unmarked_ASVbased_withUnknowns_allsamples.txt
      perl ${metapipedir}/assets/replicate_presence_absence.pl -i ${outdirectory}/ASV2Taxonomy/ASVs_counts_NOUNKNOWNS.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/presenceabsence_unmarked_ASVbased_NoUnknowns_allsamples.txt
      perl ${metapipedir}/assets/replicate_presence_absence.pl -i ${outdirectory}/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy_NOUNKNOWNS.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/presenceabsence_unmarked_TAXAbased_NoUnknowns_allsamples.txt
      perl ${metapipedir}/assets/replicate_compRelabund_detection.pl -i ${outdirectory}/dada2/ASVs_counts.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_ASVbased_withUnknowns_allsamples.txt
      perl ${metapipedir}/assets/replicate_compRelabund_detection.pl -i ${outdirectory}/ASV2Taxonomy/ASVs_counts_NOUNKNOWNS.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_ASVbased_NoUnknowns_allsamples.txt
      perl ${metapipedir}/assets/replicate_compRelabund_detection.pl -i ${outdirectory}/ASV2Taxonomy/ASVs_counts_mergedOnTaxonomy_NOUNKNOWNS.tsv -m ${outdirectory}/sample_metadata_forR.txt > ${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_TAXAbased_NoUnknowns_allsamples.txt
    fi
  fi
  mkdir ${outdirectory}/Figures/ReadsVSReplicateDetection
  if [[ "${controlspresent}" = "FALSE" && "${filterLowQualSamples}" = "FALSE" ]]; then
  Rscript --vanilla ${metapipedir}/assets/replicate_abundance_boxplot.R ${workingdirectory}/${outdirectory}/Figures/ReadsVSReplicateDetection ${workingdirectory}/${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_ASVbased_withUnknowns_allsamples.txt ${workingdirectory}/${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_ASVbased_NoUnknowns_allsamples.txt ${workingdirectory}/${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_TAXAbased_NoUnknowns_allsamples.txt \
    1>> ${workingdirectory}/${outdirectory}/Figures/ReadsVSReplicateDetection/violinboxplot_rscript_out.log 2>&1
  echo "replicate_abundance_boxplot.R	${workingdirectory}/${outdirectory}/Figures/ReadsVSReplicateDetection ${workingdirectory}/${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_ASVbased_withUnknowns_allsamples.txt ${workingdirectory}/${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_ASVbased_NoUnknowns_allsamples.txt ${workingdirectory}/${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_TAXAbased_NoUnknowns_allsamples.txt" >> ${outdirectory}/Rscript_arguments.log
  else
  Rscript --vanilla ${metapipedir}/assets/replicate_abundance_boxplot.R ${workingdirectory}/${outdirectory}/Figures/ReadsVSReplicateDetection ${workingdirectory}/${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_ASVbased_withUnknowns_filtsamples.txt ${workingdirectory}/${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_ASVbased_NoUnknowns_filtsamples.txt ${workingdirectory}/${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_TAXAbased_NoUnknowns_filtsamples.txt \
    1>> ${workingdirectory}/${outdirectory}/Figures/ReadsVSReplicateDetection/violinboxplot_rscript_out.log 2>&1
  echo "replicate_abundance_boxplot.R	${workingdirectory}/${outdirectory}/Figures/ReadsVSReplicateDetection ${workingdirectory}/${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_ASVbased_withUnknowns_filtsamples.txt ${workingdirectory}/${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_ASVbased_NoUnknowns_filtsamples.txt ${workingdirectory}/${outdirectory}/processed_tables/replicate_based_detection/compRelAbund_replicateDetection_TAXAbased_NoUnknowns_filtsamples.txt" >> ${outdirectory}/Rscript_arguments.log
  fi
  

fi #replicate if

echo "figuresFinished=TRUE" >> ${outdirectory}/progress.txt
  
fi #Final fi for if Figures folder present statement

echo "YOU MADE IT!"

sleep 1








