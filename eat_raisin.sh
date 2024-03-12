#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=batch
#SBATCH --mem-per-cpu=2gb
#SBATCH --job-name=subid

source /opt/conda/etc/profile.d/conda.sh
conda activate oatmeal_2024JAN_V2


submission=690 #711 (Lassa) # 716 (marburg) # 469 (mers) # 107 (nipah)  # can't find nipah fastqs
sub_date=2023-07-28
submitter=zcuba
ngs_ids=(
    "ILS3_2023AB05_ML-1015" # Flu
)

references=(
    "OP213369.1.fasta"

)


oatmeal_repo=/home/shared/repos/oatmeal_develop
source $oatmeal_repo/oatmeal.sh
for (( i=0; i<${#ngs_ids[@]}; i++ ));
do
    
    ngs_id="${ngs_ids[i]}"
    run_id=$(echo "${ngs_ids[i]}" | cut -c1-11)
    r1_path=/home/shared/sequencing-data/"$run_id"/"$ngs_id"_R1.fastq.gz.trimmed.fastq.gz
    r2_path=/home/shared/sequencing-data/"$run_id"/"$ngs_id"_R2.fastq.gz.trimmed.fastq.gz

    ref="${references[i]}"
    ref_no_ext="${ref%.*}"
    echo $ngs_id
    echo $ref
    output_path=/home/shared/service-data/"$submission"/"$ngs_id"_"$ref_no_ext"
	mkdir -p $output_path


    ##! Run following block if running IRMA #####################################
        
    # irma_out_path=/home/shared/service-data/"$submission"/IRMA_results/"$ngs_id"
    # mkdir -p /home/shared/service-data/"$submission"/IRMA_results
    
    # bash $oatmeal_repo/spoon.sh \
    #     -1 $r1_path \
    #     -2 $r2_path \
    #     -m "consensus" \
    #     -o "$output_path" \
    #     -n "$ngs_id" \
    #     -r /home/shared/references/"$ref" \
    #     -p "virus" \
    #     -t 16 \
    #     -i "raisin" \
    #     -y "$irma_out_path" \
    #     -k /home/shared/references/$(basename $ref .fasta).gbk

    ##!##############################################################################
    
    bash $oatmeal_repo/spoon.sh \
        -1 $r1_path \
        -2 $r2_path \
        -m "consensus" \
        -o "$output_path" \
        -n "$ngs_id" \
        -r /home/shared/references/"$ref" \
        -p "virus" \
        -t 16 \
        -i "raisin" \
        -k /home/shared/references/$(basename $ref .fasta).gbk

    SBC_Results_dir=/net/fs02.atcc.org/volume1/SBC_Data/SBC_Results/"$sub_date"_sub"$sub_id"_"$submitter"
     "mkdir -p $SBC_Results_dir/{incomplete_memos,complete_memos,complete_results}/$ngs_id"
     "chmod -R 775 $SBC_Results_dir/{incomplete_memos,complete_memos,complete_results}/$ngs_id"

     "cp $output_path/*"SbcResults"* clusterhn:"$SBC_Results_dir"/incomplete_memos/$ngs_id"
     "cp $output_path/consensus_results-notds/"$ngs_id"_reads_to_reference_"$ref_no_ext"_consensus.fasta clusterhn:"$SBC_Results_dir"/complete_results/$ngs_id"


done 
