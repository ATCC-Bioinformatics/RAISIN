#!/bin/bash
# set -e
usage() {
cat << "EOF"
    .___      .    _   _____ _ __    _
    /   \    /|    |  (      | |\   | 
    |__-'   /  \   |   `--.  | | \  | 
    |  \   /---'\  |      |  | |  \ | 
    /   \,'      \ / \___.'  / |   \| 

EOF
    echo "                                           
    Usage: 
    -o for top-level output/working directory (required),
    -f for an output name to add to the beginning of all generated files (optional, default=RAISIN_analysis),
    -t for threads (optional, default=8),

    ####### RAISIN Modes ####################
    ----------------------------------------------------------------------------------------------
    STANDARD mode:
      # Example 1: bash run_raisin.sh -m STANDARD -s SEQ -1 sample1_R1.fastq.gz -2 sample1_R2.fastq.gz -o sample1_results 
                      -f sample1 -n MN02121.1 -e username@gmail.com -d
      # Example 2: bash run_raisin.sh -m STANDARD -s VCF sample1.vcf -o sample1_results 
                      -f sample1 -r MN02121.1.fasta -g MN02121.1.gbk
      -m STANDARD,
      ****** Inputs ******
        Option 1 (if using sequencing reads):
          -s SEQ (to indicate using sequencing reads),
          -1 for filepath to forward Illumina fastq file,
          -2 for filepath to reverse Illumina fastq file

        Option 2 (if using a VCF)
          -s VCF (to indicate using a VCF)
          -v for filepath to user generated VCF

      ****** For References ******
        Option 1:
          -r for path to strain reference),
          -g for path to strain GBK,
        Option 2:
          -n for NCBI Accession ID to strain reference
          -d to download reference and gbk file based on NCBI assession ID,
          -e for Entrez email address to use to download references,
    ----------------------------------------------------------------------------------------------
    ANCHOR mode:
      # Example 1: bash run_raisin.sh -m ANCHOR -s SEQ -1 sample1_R1.fastq.gz -2 sample1_R2.fastq.gz -o sample1_results -f sample1_vs_anchor
                          -n MN02121.1 -b MZ45991.1 -e username@gmail.com -d
      # Example 2: bash run_raisin.sh -m ANCHOR -s VCF-v sample1.vcf -w anchor.vcf -o sample1_results -f sample1_vs_anchor
                          -r MN02121.1.fasta -g MN02121.1.gbk -a MZ45991.1.fasta -k MZ45991.1.gbk
      -m ANCHOR,
      ****** Inputs ******
        Option 1 (if using sequencing reads):
          -s SEQ (to indicate using sequencing reads),
          -1 for filepath to forward Illumina fastq file,
          -2 for filepath to reverse Illumina fastq file

        Option 2 (if using a VCF)
          -s VCF (to indicate using a VCF)
          -v for filepath to user generated VCF using strain reference,
          -w for filepath to user generated VCF using anchor reference

      ****** For References ******
        Option 1:
          -r for path to strain reference,
          -g for path to strain GBK,
          -a for path to anchor reference,
          -k for path to anchor GBK,
        Option 2:
          -n for NCBI Accession ID to strain reference
          -b for NCBI Accession ID to anchor reference
          -d to download reference and gbk file based on NCBI assession ID,
          -e for Entrez email address to use to download references

    ----------------------------------------------------------------------------------------------
    SERIAL_COMPARE mode:
      # Example 1: bash run_raisin.sh -m SERIAL_COMPARE -x sample1_variants.txt -y sample2_variants.txt -o comp_results -f sample1_vs_sample2
      -m standard,
      -x filepath to first variants.txt to compare,
      -y filepath to second variants.txt to compare,
    ----------------------------------------------------------------------------------------------
    
    Additional Help:
    ######## File Inputs  ####################
    If user choose either STANDARD or ANCHOR mode, there are two options for file input. Users can choose between either submitting deinterleaved Illumina
    sequencing readsets or submitting VCF files. 
    
    ######## References ####################
    If user choose either STANDARD or ANCHOR mode, a reference and a GenBank file MUST be supplied. A path to both files will suffice BUT
    users can also provide an NCBI Accession ID and an Entrez email address instead. RAISIN will then download the reference fasta and
    GenBank file based on the given NCBI Accession ID. 
    
    Option 1:
      -r/a for path to reference fasta file (required if not using -d flag),
      -g/k for path to reference gbk file (required if not using -d flag),
    Option 2:
      -n/b for NCBI Accession ID, i.e. OP213694,
      -d to download reference and gbk file based on NCBI assession ID or virus, -a or -n flags must be provided,
      -e for Entrez email address to use to download references,
    
    ################################
    and anything else for help." 1>&2
    exit 1
}

while getopts 'm:s:1:2:o:f:r:g:t:a:n:e:v:w:k:b:x:y:dh' OPTION
do
  case "$OPTION" in
    m) MODE=$OPTARG
      ;;
    s) INPUT=$OPTARG
      ;;
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
    a) ANCHOR_REF=$OPTARG
      ;;
    k) ANCHOR_GBK=$OPTARG
      ;;
    t) THREADS=$OPTARG
      ;;
    n) ACCESSION=$OPTARG
      ;;
    b) ANCHOR_ACC=$OPTARG
      ;;
    e) EMAIL=$OPTARG
      ;;
    v) VCF=$OPTARG
      ;;
    w) ANCHOR_VCF=$OPTARG
      ;;
    x) VARIANTS_PATH=$OPTARG
      ;;
    y) VARIANTS_PATH_2=$OPTARG
      ;;
    d) DOWNLOAD="true"
      ;;
    h) usage
      exit 1
      ;;
  esac
done
shift "$(($OPTIND -1))"

# #### ? Checking flags ####################################
if [ -z $MODE ] || ([ $MODE != "STANDARD" ] && [ $MODE != 'ANCHOR' ] && [ $MODE != 'SERIAL_COMPARE' ])
then
  echo "Mode (-m) must be supplied! Choose between 'STANDARD', 'ANCHOR', or 'SERIAL_COMPARE'"  
  exit
fi
if [ -z $WORKING_DIR ] # if no working directory is supplied
then
  echo "Output directory path (-o) must be supplied! Exiting pipeline..."
  exit  
fi

if [ $MODE == "STANDARD" ] || [ $MODE == "ANCHOR" ]
then
  if [ -z $INPUT ] || ([ $INPUT != "SEQ" ] && [ $INPUT != 'VCF' ]) # if no input mode was supplied
  then
    echo "Input mode (-s) must be supplied! Choose between 'SEQ' or 'VCF'"
    exit
  elif [ $INPUT == "SEQ" ]
  then
    if [ -z $FWD ] || [ -z $REV ] # if no forward/reverse reads supplied
    then
      echo "Despite being in SEQ input mode, no sequencing reads were given."
      echo "Please supply forward and reverse Illumina reads using the (-1) and (-2) flags."
      exit  
    fi
  elif [ $INPUT == "VCF" ]
  then
    if [ -z $VCF ]
    then
      echo "Despite being in VCF input mode, no VCF was given."
      echo "Please supply a VCF using the (-v) flag."
      exit  
    fi
    if [ -z $ANCHOR_VCF ] && [ $MODE == "ANCHOR" ]
    then
      echo "ANCHOR mode has been selected with VCF input mode, but no ANCHOR VCF was given."
      echo "Please supply a VCF that was generated based on the ANCHOR reference using the (-w) flag"
      exit
    fi
  fi

  if [ ! -z $DOWNLOAD ] # if users want to download a reference
  then
    if [ -z $EMAIL ] # check if they provided an email address
    then
      echo "Download flag was chosen, but no Entrez email address (-e) was given.!"
      exit
    fi
    if  [ -z $ACCESSION ] ## check if they provided an accession ID
    then
      echo "Download flag was chosen, but no NCBI accession number (-n) was given."
      exit
    fi
    if  [ -z $ANCHOR_ACCESSION ] && [ $MODE == "ANCHOR" ] ## check if they provided an accession ID for the anchor reference if in anchor mode
    then
      echo "ANCHOR mode was selected, but no NCBI accession number for the anchor (-b) was given."
      exit
    fi
  else ## if download mode was not selected
    if [ -z $REF ] && [ -z $GBK ]
    then
      echo "Since download flag was NOT chosen, a path to the reference fasta file (-r) and GBK (-g) must be provided."
      exit
    fi
    if [ -z $ANCHOR_REF ] && [ -z $ANCHOR_GBK ] && [ $MODE == "ANCHOR" ]
    then
      echo "ANCHOR mode was selected and the download flag was NOT chosen, but no ANCHOR references or gbks were given."
      echo "Please supply the path to the anchor reference fasta file (-a) and anchor GBK file (-k)"
      exit
    fi
  fi
elif [ $MODE == "SERIAL_COMPARE" ]
then
    if [ -z $VARIANTS_PATH ] && [ -z $VARIANTS_PATH_2 ]
    then
      echo "SERIAL_COMPARE mode was selected, but the paths to the two variant txt files to be compared should be provided"
      echo "Please use flags (-x) and (-y)"
      exit
    fi
fi


# # #### ? Checking flags DONE ####################################

# #### * Set main output path and make log file ####
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

# #### * Set main output path and make log file DONE ####

echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" | tee -a "$log_path" >&2
echo "$(date) RAISIN PIPELINE START!" | tee -a "$log_path" >&2
echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" | tee -a "$log_path" >&2
echo "Options passed to run_raisin.sh:" | tee -a "$log_path" >&2
echo "MODE: $MODE" | tee -a "$log_path" >&2
echo "WORKING_DIR: $WORKING_DIR" | tee -a "$log_path" >&2
echo "OUTPUT_NAME: $OUTPUT_NAME" | tee -a "$log_path" >&2
echo "THREADS: $THREADS" | tee -a "$log_path" >&2
if [ $MODE == "SERIAL_COMPARE" ]
then
  echo "VARIANT TXT 1: $VARIANTS_PATH" | tee -a "$log_path" >&2
  echo "VARIANT TXT 2: $VARIANTS_PATH_2" | tee -a "$log_path" >&2
else
  echo "INPUT MODE: $INPUT" | tee -a "$log_path" >&2
  echo "FWD: $FWD" | tee -a "$log_path" >&2
  echo "REV: $REV" | tee -a "$log_path" >&2
  echo "STRAIN VCF: $VCF" | tee -a "$log_path" >&2
  echo "REFERENCE PATH: $REF" | tee -a "$log_path" >&2
  echo "GBK PATH: $GBK" | tee -a "$log_path" >&2
  echo "ANCHOR REFERENCE PATH: $ANCHOR_REF" | tee -a "$log_path" >&2
  echo "ANCHOR GBK PATH: $ANCHOR_GBK" | tee -a "$log_path" >&2
  echo "DOWNLOAD REFERENCES: $DOWNLOAD" | tee -a "$log_path" >&2
  echo "ACCESSION: $ACCESSION" | tee -a "$log_path" >&2
  echo "ANCHOR ACCESSION: $ANCHOR_ACC" | tee -a "$log_path" >&2
  echo "ORGANISM NAME: $ORG_NAME" | tee -a "$log_path" >&2
  echo "EMAIL: $EMAIL" | tee -a "$log_path" >&2
fi
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
# Set output name of fasta and gbk based on user inputs
download_name="$ACCESSION"

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

#### ! Check if all file paths exist is present ###########################

if [ $MODE == "SERIAL_COMPARE" ]
then
  file_check=($VARIANTS_PATH $VARIANTS_PATH_2)
elif [ $MODE == "STANDARD" ]
then
  if [ $INPUT == "SEQ" ]
  then
    file_check=($FWD $REV $REF $GBK)
  else
    file_check=($VCF $REF $GBK)
  fi
elif [ $MODE == "ANCHOR" ]
then
  if [ $INPUT == "SEQ" ]
  then
    file_check=($FWD $REV $REF $GBK $ANCHOR_REF $ANCHOR_GBK)
  else
    file_check=($VCF $ANCHOR_VCF $REF $GBK $ANCHOR_REF $ANCHOR_GBK)
  fi
fi

for (( i=0; i<${#file_check[@]}; i++ ));
do
    file_to_check="${file_check[i]}"
    if [ ! -f $file_to_check ]
    then
      echo "$(date) $file_to_check does not exist! Please verify file path" | tee -a "$log_path" >&2
      exit
    fi
done



# #### ! Reference Download DONE ####################################
# runID=$(echo $FWD | cut -c1-11)

# echo "Running consensus pipeline" | tee -a "$log_path" >&2
# consensus \
#   $WORKING_DIR \
#   $tag \
#   $NGS_ID \
#   $FWD \
#   $REV \
#   $THREADS \
#   $nucID \
#   $REF


# ref_id=$(basename $REF .fasta)
# outdir=$output_path/consensus_results"$tag"
# vcf_path=$outdir/"$NGS_ID"_"reads_to_reference"_"$ref_id"_lofreq.vcf
# variants_txt=$outdir/"$NGS_ID"_"reads_to_reference"_"$ref_id"_lofreq_variants.txt

# if [ ! -f $variants_txt ]
# then
#   if [ ! -f $vcf_path ]
#   then
#     echo "$(date) VCF does not exist, check consensus results logs" | tee -a "$log_path" >&2
#     continue
#   else
#     logger "creation of variants.txt for raisin" "Start" | tee -a "$log_path" >&2 
#     source /opt/conda/etc/profile.d/conda.sh
#     conda init bash
#     conda deactivate
#     ### activate raisin environment and run raisin
#     source /home/shared/service-data/MERS/feature_viewer_env/bin/activate
#     python $oatmeal_path/scripts_for_raisin/universal_raisin.py "$vcf_path" $REF $GBK
#     ### deactivate raisin environment
#     deactivate
#     source /opt/conda/etc/profile.d/conda.sh
#     conda activate oatmeal_2024JAN_V2.bak

#     if [ -f $variants_txt ]
#     then
#       logger "creation of variants.txt for raisin" "True" | tee -a "$log_path" >&2 
#     else
#       logger "creation of variants.txt for raisin" "Fail" | tee -a "$log_path" >&2 
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
#   #! ANCHOR mode must providing sequencing reads
#   #! anchor mode: map reads to strain reference, call variants, call consensus, mafft, retrive variants.txt
#   reads_to_reference_consensus=""
#   reads_to_reference_final_vcf=""

#   cat $REF > $WORKING_DIR/"$OUTPUT_NAME"_mafft_input.fasta
#   cat $ANCHOR_REF >> $WORKING_DIR/"$OUTPUT_NAME"_mafft_input.fasta
#   cat $reads_to_reference_consensus >> $WORKING_DIR/"$OUTPUT_NAME"_mafft_input.fasta
#   mafft $WORKING_DIR/"$OUTPUT_NAME"_mafft_input.fasta > $WORKING_DIR/"$OUTPUT_NAME"_mafft_alignment.fasta

#   # Iterate down alignment and determine variant type, genetic region, codon mutation, etc.
#   variant_script_path=$(dirname "$0")/BEI_variants.py
#   python3 $variant_script_path -f $reference_dir_path/"$OUTPUT_NAME"_mafft_alignment.fasta -v $reads_to_reference_final_vcf -o $WORKING_DIR/variants.txt