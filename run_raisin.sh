#!/bin/bash
# set -e
usage() {
    echo "Usage: run appropriate pipeline for submission
    -1 for forward Illumina fastq file (required),
    -2 for reverse Illumina fastq file (required),
    -o for top-level output/working directory (required),
    -f for an output name to add to the beginning of all generated files,
    -r for path to reference fasta file (required if not using -d flag),
    -g for path to reference gbk file (required if not using -d flag),
    -a for NCBI Accession ID, i.e. OP213694,
    -n for virus name,
    -d to download reference and gbk file based on NCBI assession ID or virus, -a or -n flags must be provided,
    -e for Entrez email address to use to download references,
    -t for threads,
    and anything else for help." 1>&2
    exit 1
}

while getopts '1:2:o:f:r:g:t:a:n:e:d' OPTION
do
  case "$OPTION" in
    1) FWD=$OPTARG
      ;;
    2) REV=$OPTARG
      ;;
    o) WORKING_DIR=$OPTARG
      ;;
    f) OUTPUT_NAME=$OPTARG
      ;;
    r) REF=$OPTARG
      ;;
    g) GBK=$OPTARG
      ;;
    t) THREADS=$OPTARG
      ;;
    a) ACCESSION=$OPTARG
      ;;
    n) ORG_NAME=$OPTARG
      ;;
    e) EMAIL=$OPTARG
      ;;
    d) DOWNLOAD="true"
      ;;
    ?) usage
      exit 1
      ;;
  esac
done
shift "$(($OPTIND -1))"

#### ? Checking flags ####################################
if [ -z $FWD ] || [ -z $REV ] # if no forward/reverse reads supplied
then
  echo "Forward and reverse Illumina reads must be supplied! Exiting pipeline"
  exit  
fi
if [ -z $WORKING_DIR ] # if no working directory is supplied
then
  echo "Output directory path must be supplied! Exiting pipeline"
  exit  
fi
if [ -z $REF ] || [ -z $GBK ] # if a reference or gbk is not supllied
then
    if [ ! -z $DOWNLOAD ] # if users want to download a reference
    then
      if [ -z $EMAIL ]
      then
        echo "An Entrez email address must be given in order to download sequences from NCBI! Exiting pipeline"
        exit
      fi
      if  [ ! -z $ACCESSION ] && [ ! -z $ORG_NAME ] # check if user supplied both an accession number + organism name
      then
        echo "Choose one! If download flag is chosen, choose either NCBI accession number (-a) or an organism name (-n)! Exiting pipeline"
        exit
      elif [ -z $ACCESSION ] && [ -z $ORG_NAME ] # check if no accession number or organism name was supplied
      then
        echo "If download flag is chosen, an NCBI accession number (-a) or an organism name (-n) must be given! Exiting pipeline"
        exit
      fi
    fi 
    if [ -z $DOWNLOAD ] && ([ -z $REF ] || [ -z $GBK ]) # if user did not want to download a reference but did not supply a ref/gbk path
    then
      echo "A reference and gbk path must be given if not using the download flag! Exiting pipeline"
      exit
    fi 
fi

if [ -z $OUTPUT_NAME ] # if no output name is given, use a default
then
  OUTPUT_NAME="RAISIN_analysis"
fi

#### ? Checking flags DONE ####################################

#### * Set main output path and make log file ####
output_path="$WORKING_DIR"
log_path="$WORKING_DIR"/"$OUTPUT_NAME".log
if [ ! -d "$output_path" ]
then
  mkdir -p $output_path
  chmod -R 775 "$output_path"
fi
if [ ! -f "$log_path" ]
then
  touch "$log_path"
  chmod 775 $log_path
fi

### REDIRECT STDOUT AND STDERR ####
exec 1> "$output_path"/"$(date +"%Y_%b_%d")"_"$OUTPUT_NAME"_stdOUT.log
chmod 775 "$output_path"/"$(date +"%Y_%b_%d")"_"$OUTPUT_NAME"_stdOUT.log
exec 2> "$output_path"/"$(date +"%Y_%b_%d")"_"$OUTPUT_NAME"_stdERR.log
chmod 775 "$output_path"/"$(date +"%Y_%b_%d")"_"$OUTPUT_NAME"_stdERR.log

#### * Set main output path and make log file DONE ####

echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" | tee -a "$log_path" >&2
echo "$(date) RAISIN PIPELINE START!" | tee -a "$log_path" >&2
echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" | tee -a "$log_path" >&2
echo "Options passed to run_raisin.sh:" | tee -a "$log_path" >&2
echo "FWD: $FWD" | tee -a "$log_path" >&2
echo "REV: $REV" | tee -a "$log_path" >&2
echo "WORKING_DIR: $WORKING_DIR" | tee -a "$log_path" >&2
echo "OUTPUT_NAME: $OUTPUT_NAME" | tee -a "$log_path" >&2
echo "REFERENCE PATH: $REF" | tee -a "$log_path" >&2
echo "GBK PATH: $GBK" | tee -a "$log_path" >&2
echo "THREADS: $THREADS" | tee -a "$log_path" >&2
echo "DOWNLOAD REFERENCES: $DOWNLOAD" | tee -a "$log_path" >&2
echo "ACCESSION: $ACCESSION" | tee -a "$log_path" >&2
echo "ORGANISM NAME: $ORG_NAME" | tee -a "$log_path" >&2
echo "EMAIL: $EMAIL" | tee -a "$log_path" >&2
echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" | tee -a "$log_path" >&2

logger(){
    #written for general cleanliness in other functions
    phrase=$1
    action_status=$2
    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
    if [ $action_status = "Success" ]; then
        echo "$(date) $phrase was successful!"| tee -a "$log_path" >&2
    elif [ $action_status == "Start" ]; then
        echo "$(date) Starting $phrase "| tee -a "$log_path" >&2
    elif [ $action_status == "Done" ]; then
        echo "$(date) Y'all've already done this. Passing $phrase"| tee -a "$log_path" >&2
    elif [ $action_status == "Fail" ]; then
        echo "$(date) $phrase failed."| tee -a "$log_path" >&2
    else
        echo "$(date) You forgot a status - this message should never be seen"| tee -a "$log_path" >&2
    fi
    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
}
export -f logger

#### ! Reference Download ####################################
## Set output name of fasta and gbk based on user inputs
if [ ! -z $ACCESSION ] # if user requesed download via accession
then
  download_name="$ACCESSION"
elif [ ! -z $ORG_NAME ] ## if user requesed download via organism name
then
  download_name="$ORG_NAME"
fi


if [ ! -z $DOWNLOAD ] # if user wants to download files from NCBI
then
  REF="$output_path"/"$download_name".fasta
  GBK="$output_path"/"$download_name".gbk
  if [ ! -f $REF ] && [ ! -f $GBK ] # if fasta and gbk haven't been downloaded
  then
    logger "reference and GBK file download from NCBI using $download_name" "Start" | tee -a "$log_path" >&2
    if [ ! -z $ACCESSION ] # if user requesed download via accession
    then
      python /home/shared/repos/RAISIN_develop/get_accs.py $ACCESSION $EMAIL $output_path $log_path
    elif [ ! -z $ORG_NAME ] ## if user requesed download via organism name
    then
      echo "Need to code reference download based on organism name"
      # grab organism name from assembly_summary.txt? Grab representative genome
      # make sure fasta and gbk download is not 0kb
    fi

    if [ ! -f $REF ] && [ ! -f $GBK ]
    then
      logger "reference and GBK file download from NCBI using $download_name" "Fail" | tee -a "$log_path" >&2
      echo "$(date) Make sure $ACCESSION and $EMAIL or $ORG_NAME are valid." | tee -a "$log_path" >&2
    else
      logger "reference and GBK file download from NCBI using $download_name" "Success" | tee -a "$log_path" >&2
    fi
  else
    logger "reference and GBK file download from NCBI using $download_name" "DONE" | tee -a "$log_path" >&2
  fi
fi

#### ! Reference Download DONE ####################################

# output_path=$1
# tag=$2
# NGS_ID=$3
# PIPELINE=$4
# FWD=$5
# REV=$6
# THREADS=$7
# nucID=$8
# main_commit=$9
# MODE=${10}
# REF=${11}
# oatmeal_path=${12}
# GBK=${13}
# TEST=${14}
# DEBUG=${15}

# runID=$(echo $FWD | cut -c1-11)

# echo "Running consensus pipeline" | tee -a "$log_path" >&2
#   if [ -z $REF ]
#   then
#     echo "$(date) Consensus called without providing a reference " | tee -a "$log_path" >&2
#     continue
#   else
#     consensus \
#       $WORKING_DIR \
#       $tag \
#       $NGS_ID \
#       $PIPELINE \
#       $FWD \
#       $REV \
#       $THREADS \
#       $nucID \
#       $main_commit \
#       $MODE \
#       $REF \
#       $TEST \
#       $DEBUG
#   fi
# else
#     echo "Please select an analysis mode from de_novo or consensus" | tee -a "$log_path" >&2
# fi

# ref_id=$(basename $REF .fasta)
# outdir=$output_path/consensus_results"$tag"
# vcf_path=$outdir/"$NGS_ID"_"reads_to_reference"_"$ref_id"_lofreq.vcf
# variants_txt=$outdir/"$NGS_ID"_"reads_to_reference"_"$ref_id"_lofreq_variants.txt
# if [ ! -f $variants_txt ]
# then
#   if [ ! -f $GBK ]
#   then
#     echo "$(date) $GBK does not exist" | tee -a "$log_path" >&2
#     continue
#   else
#     if [ ! -f $vcf_path ]
#     then
#       echo "$(date) VCF does not exist, check consensus results logs" | tee -a "$log_path" >&2
#       continue
#     else
#       logger "creation of variants.txt for raisin" "Start" | tee -a "$log_path" >&2 
#       source /opt/conda/etc/profile.d/conda.sh
#       conda init bash
#       conda deactivate
#       ### activate raisin environment and run raisin
#       source /home/shared/service-data/MERS/feature_viewer_env/bin/activate
#       python $oatmeal_path/scripts_for_raisin/universal_raisin.py "$vcf_path" $REF $GBK
#       ### deactivate raisin environment
#       deactivate
#       source /opt/conda/etc/profile.d/conda.sh
#       conda activate oatmeal_2024JAN_V2.bak

#       if [ -f $variants_txt ]
#       then
#         logger "creation of variants.txt for raisin" "True" | tee -a "$log_path" >&2 
#       else
#         logger "creation of variants.txt for raisin" "Fail" | tee -a "$log_path" >&2 
#       fi
#     fi
#   fi
# else
#   logger "creation of variants.txt for raisin" "Pass" | tee -a "$log_path" >&2 
# fi

# docx=$( ls "$output_path"/*.docx 2>/dev/null ) # 2>/dev/null mutes the error message if the file can't be found
# if [ -z $docx ]
# then
#     echo "$(date) Creating report in 'consensus' mode." | tee -a "$log_path" >&2
#     /home/src/anaconda3/envs/oatmeal_2024JAN_V2/bin/Rscript $oatmeal_path/bowl.R -o $WORKING_DIR -n $NGS_ID -r $REF -c -p $PIPELINE -b $oatmeal_path -v $variants_txt
# else
#   echo "$(date) Report already created" | tee -a "$log_path" >&2
# fi
# #     echo "$(date) Report was not created, maybe a sample.json issue. Manually creating it" | tee -a "$log_path" >&2
# #     python $oatmeal_path/scripts_for_raisin/create_sample_json.py \
# #         $oatmeal_path/scripts_for_raisin/sample_template.json \
# #         "$output_path"/pre_QC/fastp.json \
# #         "$NGS_ID" \
# #         "$output_path"

# #     Rscript $oatmeal_path/bowl.R -o $WORKING_DIR -n $NGS_ID -r $REF -c -p $PIPELINE -b $oatmeal_path
# # fi