#!/bin/bash
# set -e
shopt -s expand_aliases
source /home/shared/bash.bashrc
source /opt/conda/etc/profile.d/conda.sh

logger(){
    #written for general cleanliness in other functions
    cmd=$1
    pass=$2
    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
    if [ $pass = "True" ]; then
        echo "$(date) $cmd was successful!"| tee -a "$log_path" >&2
    elif [ $pass == "Start" ]; then
        echo "$(date) Starting $cmd "| tee -a "$log_path" >&2
    elif [ $pass == "Pass" ]; then
        echo "$(date) Y'all've already done this. Passing $cmd"| tee -a "$log_path" >&2
    elif [ $pass == "Fail" ]; then
        echo "$(date) $cmd failed. Please check error log for details."| tee -a "$log_path" >&2
    else
        echo "$(date) You forgot a status - this message should never be seen"| tee -a "$log_path" >&2
    fi
    echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------" >&2
}
export -f logger


run_raisin(){
  output_path=$1
  tag=$2
  NGS_ID=$3
  PIPELINE=$4
  FWD=$5
  REV=$6
  THREADS=$7
  nucID=$8
  main_commit=$9
  MODE=${10}
  REF=${11}
  oatmeal_path=${12}
  GBK=${13}
  TEST=${14}
  DEBUG=${15}

  runID=$(echo $FWD | cut -c1-11)

  echo "MODE: $MODE" | tee -a "$log_path" >&2
  echo "REF: $REF" | tee -a "$log_path" >&2
  echo "Oatmeal path: $oatmeal_path" | tee -a "$log_path" >&2
  echo "GBK: $GBK" | tee -a "$log_path" >&2
  if [ "$MODE" == "consensus" ]
  then
    # tag="-notds"
    echo "Running consensus pipeline" | tee -a "$log_path" >&2
    if [ -z $REF ]
    then
      echo "$(date) Consensus called without providing a reference " | tee -a "$log_path" >&2
      continue
    else
      consensus \
        $WORKING_DIR \
        $tag \
        $NGS_ID \
        $PIPELINE \
        $FWD \
        $REV \
        $THREADS \
        $nucID \
        $main_commit \
        $MODE \
        $REF \
        $TEST \
        $DEBUG
    fi
  else
      echo "Please select an analysis mode from de_novo or consensus" | tee -a "$log_path" >&2
  fi

  ref_id=$(basename $REF .fasta)
  outdir=$output_path/consensus_results"$tag"
  vcf_path=$outdir/"$NGS_ID"_"reads_to_reference"_"$ref_id"_lofreq.vcf
  variants_txt=$outdir/"$NGS_ID"_"reads_to_reference"_"$ref_id"_lofreq_variants.txt
  if [ ! -f $variants_txt ]
  then
    if [ ! -f $GBK ]
    then
      echo "$(date) $GBK does not exist" | tee -a "$log_path" >&2
      continue
    else
      if [ ! -f $vcf_path ]
      then
        echo "$(date) VCF does not exist, check consensus results logs" | tee -a "$log_path" >&2
        continue
      else
        logger "creation of variants.txt for raisin" "Start" | tee -a "$log_path" >&2 
        source /opt/conda/etc/profile.d/conda.sh
        conda init bash
        conda deactivate
        ### activate raisin environment and run raisin
        source /home/shared/service-data/MERS/feature_viewer_env/bin/activate
        python $oatmeal_path/scripts_for_raisin/universal_raisin.py "$vcf_path" $REF $GBK
        ### deactivate raisin environment
        deactivate
        source /opt/conda/etc/profile.d/conda.sh
        conda activate oatmeal_2024JAN_V2.bak

        if [ -f $variants_txt ]
        then
          logger "creation of variants.txt for raisin" "True" | tee -a "$log_path" >&2 
        else
          logger "creation of variants.txt for raisin" "Fail" | tee -a "$log_path" >&2 
        fi
      fi
    fi
  else
    logger "creation of variants.txt for raisin" "Pass" | tee -a "$log_path" >&2 
  fi

  docx=$( ls "$output_path"/*.docx 2>/dev/null ) # 2>/dev/null mutes the error message if the file can't be found
  if [ -z $docx ]
  then
      echo "$(date) Creating report in 'consensus' mode." | tee -a "$log_path" >&2
      /home/src/anaconda3/envs/oatmeal_2024JAN_V2/bin/Rscript $oatmeal_path/bowl.R -o $WORKING_DIR -n $NGS_ID -r $REF -c -p $PIPELINE -b $oatmeal_path -v $variants_txt
  else
    echo "$(date) Report already created" | tee -a "$log_path" >&2
  fi
  #     echo "$(date) Report was not created, maybe a sample.json issue. Manually creating it" | tee -a "$log_path" >&2
  #     python $oatmeal_path/scripts_for_raisin/create_sample_json.py \
  #         $oatmeal_path/scripts_for_raisin/sample_template.json \
  #         "$output_path"/pre_QC/fastp.json \
  #         "$NGS_ID" \
  #         "$output_path"

  #     Rscript $oatmeal_path/bowl.R -o $WORKING_DIR -n $NGS_ID -r $REF -c -p $PIPELINE -b $oatmeal_path
  # fi
}