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
    -q for minimum variant frequency threshold, enter as a decimal, i.e. 0.05, (optional, default=0.05),
    -c for minimum read coverage for variant calling, (optional, default=10) 

    ####### RAISIN Modes ####################
    ----------------------------------------------------------------------------------------------
    STANDARD mode:
      # Example 1: bash raisin.sh -m STANDARD -s SEQ -1 sample1_R1.fastq.gz -2 sample1_R2.fastq.gz -o sample1_results 
                      -f sample1 -n MN02121.1 -e username@gmail.com -d
      # Example 2: bash raisin.sh -m STANDARD -s VCF -v sample1.vcf -o sample1_results 
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
      # Example 1: bash raisin.sh -m ANCHOR -s SEQ -1 sample1_R1.fastq.gz -2 sample1_R2.fastq.gz -o sample1_results -f sample1_vs_anchor
                          -n MN02121.1 -b MZ45991.1 -e username@gmail.com -d
      # Example 2: bash raisin.sh -m ANCHOR -s SEQ -1 sample1_R1.fastq.gz -2 sample1_R2.fastq.gz -o sample1_results -f sample1_vs_anchor
                          -r MN02121.1.fasta -a MZ45991.1.fasta -k MZ45991.1.gbk
      -m ANCHOR,
      ****** Inputs ******
          -s SEQ (to indicate using sequencing reads),
          -1 for filepath to forward Illumina fastq file,
          -2 for filepath to reverse Illumina fastq file

      ****** For References ******
        Option 1:
          -r for path to strain reference,
          -a for path to anchor reference,
          -k for path to anchor GBK,
        Option 2:
          -n for NCBI Accession ID to strain reference
          -b for NCBI Accession ID to anchor reference
          -d to download reference and gbk file based on NCBI assession ID,
          -e for Entrez email address to use to download references

    ----------------------------------------------------------------------------------------------
    COMPARE mode:
      # Example 1: bash raisin.sh -m COMPARE -x sample1_variants.txt -y sample2_variants.txt 
                    -o comp_results -f sample1_vs_sample2 -X Sample_1 -Y Sample_2
      -m COMPARE,
      -x filepath to first variants.txt to compare,
      -y filepath to second variants.txt to compare,
      -X identifier for first variants.txt, i.e. Passage_0
      -Y identifier for second variants.txt, i.e. Passage_1
    ----------------------------------------------------------------------------------------------
    
    Additional Help:
    ######## File Inputs  ####################
    If the user chooses STANDARD mode, there are two options for file input. Users can choose between either submitting deinterleaved Illumina
    sequencing readsets or submitting VCF files. 
    
    ######## References ####################
    If user choose either STANDARD or ANCHOR mode, a reference and a GenBank file MUST be supplied. A path to both files will suffice BUT
    users can also provide an NCBI Accession ID and an Entrez email address instead. RAISIN will then download the reference fasta and
    GenBank file based on the given NCBI Accession ID. 
    
    Option 1:
      -r/a for path to reference fasta file (required if not using -d flag), (use -a for reference path for anchor)
      -g/k for path to reference gbk file (required if not using -d flag), (use -k for annotations path for anchor)
    Option 2:
      -n/b for NCBI Accession ID, i.e. OP213694,
      -d to download reference and gbk file based on NCBI assession ID or virus, -a or -n flags must be provided,
      -e for Entrez email address to use to download references,
    
    ################################
    and anything else for help." 1>&2
    exit 1
}

while getopts 'm:s:1:2:o:f:r:g:t:q:c:a:n:e:v:w:k:b:x:y:X:Y:dh' OPTION
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
    q) MIN_FREQ=$OPTARG
      ;;
    c) MIN_COV=$OPTARG
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
    X) SAMPLE_NAME_1=$OPTARG
      ;;
    Y) SAMPLE_NAME_2=$OPTARG
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
if [ -z $MODE ] || ([ $MODE != "STANDARD" ] && [ $MODE != 'ANCHOR' ] && [ $MODE != 'COMPARE' ])
then
  echo "Mode (-m) must be supplied! Choose between 'STANDARD', 'ANCHOR', or 'COMPARE'"  
  exit
fi
if [ -z $WORKING_DIR ] # if no working directory is supplied
then
  echo "Output directory path (-o) must be supplied! Exiting pipeline..."
  exit  
fi
## Setting default values if necessary
if [ -z $THREADS ] # if no threads is supplied
then
  echo "Setting threads to default 8"
  THREADS=8
fi
if [ -z $MIN_FREQ ] # if no minium variant frequency is supplied
then
  echo "Setting minimum frequency to default 0.05"
  MIN_FREQ=0.05
fi
if [ -z $MIN_COV ] # if no minimum read coverage is supplied
then
  echo "Setting minimum coverage to default 10"
  MIN_COV=10
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
  elif [ $INPUT == "VCF" ] && [ $MODE == "STANDARD" ]
  then
    if [ -z $VCF ]
    then
      echo "Despite being in VCF input mode, no VCF was given."
      echo "Please supply a VCF using the (-v) flag."
      exit  
    fi
  elif [ $INPUT == "VCF" ] && [ $MODE == "ANCHOR" ]
  then
    echo "ANCHOR mode cannot be run in VCP input mode, please only run in SEQ input mode"
    exit
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
    if [ -z $REF ] && [ -z $GBK ] && [ $MODE == "STANDARD" ]
    then
      echo "Since download flag was NOT chosen, a path to the reference fasta file (-r) and GBK (-g) must be provided for STANDARD mode."
      exit
    fi
    if [ -z $REF ] && [ $MODE == "ANCHOR" ]
    then
      echo "Since download flag was NOT chosen, a path to the reference fasta file (-r) must be provided for ANCHOR mode."
      exit
    fi
    if [ -z $ANCHOR_REF ] && [ -z $ANCHOR_GBK ] && [ $MODE == "ANCHOR" ]
    then
      echo "ANCHOR mode was selected and the download flag was NOT chosen, but no ANCHOR references or gbks were given."
      echo "Please supply the path to the anchor reference fasta file (-a) and anchor GBK file (-k)"
      exit
    fi
  fi
elif [ $MODE == "COMPARE" ]
then
    if [ -z $VARIANTS_PATH ] && [ -z $VARIANTS_PATH_2 ]
    then
      echo "COMPARE mode was selected, but the paths to the two variant txt files to be compared should be provided"
      echo "Please use flags (-x) and (-y)"
      exit
    fi
fi


# # #### ? Checking flags DONE ####################################

# #### * Set main output path and make log file ####
output_path="$WORKING_DIR"
log_path="$WORKING_DIR"/"$OUTPUT_NAME".log
outdir="$WORKING_DIR"/raisin_results_"$MODE"
source $(dirname $0)/raisin_functions.sh

if [ ! -d "$output_path" ]
then
  mkdir -p $output_path
  chmod -R 775 "$output_path"
fi
if [ ! -d "$outdir" ]
then
  mkdir -p $outdir
  chmod -R 775 "$outdir"
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
if [ $MODE == "COMPARE" ]
then
  echo "VARIANT TXT 1: $VARIANTS_PATH" | tee -a "$log_path" >&2
  echo "VARIANT TXT 2: $VARIANTS_PATH_2" | tee -a "$log_path" >&2
  echo "VARIANT 1 IDENTIFIER: $SAMPLE_NAME_1" | tee -a "$log_path" >&2
  echo "VARIANT 2 IDENTIFIER: $SAMPLE_NAME_2" | tee -a "$log_path" >&2
else
  echo "INPUT MODE: $INPUT" | tee -a "$log_path" >&2
  echo "FWD: $FWD" | tee -a "$log_path" >&2
  echo "REV: $REV" | tee -a "$log_path" >&2
  echo "MINIMUM VARIANT FREQUENCY: $MIN_FREQ" | tee -a "$log_path" >&2
  echo "MINIMUM READ COVERAGE: $MIN_COV" | tee -a "$log_path" >&2
  echo "STRAIN VCF (STANDARD ONLY): $VCF" | tee -a "$log_path" >&2
  echo "REFERENCE PATH: $REF" | tee -a "$log_path" >&2
  echo "GBK PATH: $GBK" | tee -a "$log_path" >&2
  echo "DOWNLOAD REFERENCES: $DOWNLOAD" | tee -a "$log_path" >&2
  echo "ACCESSION: $ACCESSION" | tee -a "$log_path" >&2
  echo "ORGANISM NAME: $ORG_NAME" | tee -a "$log_path" >&2
  echo "EMAIL: $EMAIL" | tee -a "$log_path" >&2
  if [ $MODE == "ANCHOR" ]
  then
    echo "ANCHOR REFERENCE PATH: $ANCHOR_REF" | tee -a "$log_path" >&2
    echo "ANCHOR GBK PATH: $ANCHOR_GBK" | tee -a "$log_path" >&2
    echo "ANCHOR ACCESSION: $ANCHOR_ACC" | tee -a "$log_path" >&2
  fi
fi
echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" | tee -a "$log_path" >&2


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
      python $(dirname $0)/get_accs.py $ACCESSION $EMAIL $output_path $log_path
    fi

    if [ ! -f $REF ] && [ ! -f $GBK ]
    then
      logger "reference and GBK file download from NCBI using $download_name" "Fail" | tee -a "$log_path" >&2
      echo "$(date) Make sure $ACCESSION and $EMAIL or $ORG_NAME are valid." | tee -a "$log_path" >&2
    else
      logger "reference and GBK file download from NCBI using $download_name" "Success" | tee -a "$log_path" >&2
    fi
  else
    logger "reference and GBK file download from NCBI using $download_name" "Done" | tee -a "$log_path" >&2
  fi
fi
# #### ! Reference Download DONE ####################################

#### ! Check if all file paths exist is present ###########################

if [ $MODE == "COMPARE" ]
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
    file_check=($FWD $REV $REF $ANCHOR_REF $ANCHOR_GBK)
  else
    file_check=($VCF $REF $ANCHOR_REF $ANCHOR_GBK)
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


echo "$(date) Running in $MODE mode..." | tee -a "$log_path" >&2
if [ $MODE == "COMPARE" ]
then
  python $(dirname $0)/combine_variant_txt.py $VARIANTS_PATH $VARIANTS_PATH_2 $SAMPLE_NAME_1 $SAMPLE_NAME_2 $outdir
else
  if [ $INPUT == "SEQ" ]
  then
    # trim short reads
    trim_filt "$WORKING_DIR" $THREADS $FWD $REV

    if (file "$FWD" | grep -q compressed ) ; then  # if file is compressed, then use different basename extension
      FWD="$WORKING_DIR"/$(basename $FWD .fastq.gz).filtered.fastq.gz
      REV="$WORKING_DIR"/$(basename $REV .fastq.gz).filtered.fastq.gz
    else
        FWD="$WORKING_DIR"/$(basename $FWD .fastq).filtered.fastq.gz
        REV="$WORKING_DIR"/$(basename $REV .fastq).filtered.fastq.gz
    fi
    # run multiqc on trimmed reads
    run_multiqc "$WORKING_DIR" "$FWD" "$REV" "$THREADS"
    consensus \
        $outdir \
        $OUTPUT_NAME \
        $FWD \
        $REV \
        $THREADS \
        $REF \
        $MIN_FREQ \
        $MIN_COV
    
    ref_id=$(basename $REF .fasta)
    vcf_path=$outdir/"$OUTPUT_NAME"_"reads_to_reference"_"$ref_id"_lofreq.vcf
    variants_txt=$outdir/"$OUTPUT_NAME"_"reads_to_reference"_"$ref_id"_lofreq_variants.txt

    if [ $MODE == "ANCHOR" ] #! Does ANCHOR only work with reads?? (Yes!)
    then
      reads_to_reference_consensus=$outdir/"$OUTPUT_NAME"_"reads_to_reference"_"$ref_id"_consensus.fasta
      cat $ANCHOR_REF > $outdir/"$OUTPUT_NAME"_mafft_input.fasta
      cat $REF >> $outdir/"$OUTPUT_NAME"_mafft_input.fasta
      cat $reads_to_reference_consensus >> $outdir/"$OUTPUT_NAME"_mafft_input.fasta
      mafft $outdir/"$OUTPUT_NAME"_mafft_input.fasta > $outdir/"$OUTPUT_NAME"_mafft_alignment.fasta
      
      mafft_file=$outdir/"$OUTPUT_NAME"_mafft_alignment.fasta
      python $(dirname $0)/universal_raisin.py \
        $vcf_path \
        $ANCHOR_GBK \
        $MODE \
        $OUTPUT_NAME \
        $outdir \
        $mafft_file
    
    elif [ $MODE == "STANDARD" ]
    then
      python $(dirname $0)/universal_raisin.py \
        $vcf_path \
        $GBK \
        $MODE \
        $OUTPUT_NAME \
        $outdir
    fi
  
  elif [ $INPUT == "VCF" ] && [ $MODE == "STANDARD" ]
  then
    python $(dirname $0)/universal_raisin.py \
        $VCF \
        $GBK \
        $MODE \
        $OUTPUT_NAME \
        $outdir
  fi
fi


echo "$(date) RAISIN completed successfully!" | tee -a "$log_path" >&2
