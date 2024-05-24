#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=normal
#SBATCH --job-name=raisin_test

source /opt/conda/etc/profile.d/conda.sh
conda activate oatmeal_2024JAN_V2

ngs_ids=(
    "ILS3_2023AB05_ML-1015" # Standard
)

references=(
    "OP213369.1.fasta" # Reference download 
)


# vcf="/home/shared/service-data/863/ILM3_2023CF07_NR-59685_EPI_ISL_18507690/reads_to_reference/ILM3_2023CF07_NR-59685_reads_to_reference_EPI_ISL_18507690_final.vcf"
# mafft="/home/shared/service-data/863/ILM3_2023CF07_NR-59685_EPI_ISL_18507690/reads_to_reference/ILM3_2023CF07_NR-59685_mafft_alignment.fasta"
# gbk="/home/shared/references/NC_045512.2.gbk"

# # python msa_raisin.py
# python universal_raisin.py $vcf $gbk "ANCHOR" "full_test" $outdir $mafft

working_dir="/home/shared/projects/RAISIN/testing"
#########? STANDARD TEST: SEQ MODE ######################## (with reference download)
ngs_ids=(
    "ILS3_2023AB05_ML-1015" # Standard
    "ILM3_2024AH04_229E-HTB-55_P3" # submission 955
)

references=(
    "OP213369.1.fasta"
    "MW202340.1.fasta"
)


# for (( i=0; i<${#ngs_ids[@]}; i++ ));
# do
    
#     ngs_id="${ngs_ids[i]}"
#     run_id=$(echo "${ngs_ids[i]}" | cut -c1-11)
#     r1_path=/home/shared/raw_fastqs/"$run_id"/"$ngs_id"_R1.fastq.gz
#     r2_path=/home/shared/raw_fastqs/"$run_id"/"$ngs_id"_R2.fastq.gz

#     ref="${references[i]}"
#     ref_no_ext="${ref%.*}" 
#     output_dir="$working_dir"/STANDARD/"$ngs_id"_seq_mode
#     bash /home/shared/repos/RAISIN_develop/run_raisin.sh \
#         -m "STANDARD" \
#         -s "SEQ" \
#         -1 $r1_path \
#         -2 $r2_path \
#         -o $output_dir \
#         -f "$ngs_id" \
#         -t 16 \
#         -n "$ref_no_ext" \
#         -e "dyarmosh@atcc.org" \
#         -d

# done 

######################################################?

#########? STANDARD TEST: VCF MODE ######################## (with NO reference download)
vcfs=(
    "/home/shared/service-data/994/ILM3_2024AN04_NR-540_AF528167/consensus_results-notds/ILM3_2024AN04_NR-540_reads_to_reference_AF528167_lofreq.vcf"
    "/home/shared/service-data/955/ILM3_2024AH06_229E-HTB-46_P3_MW202340.1/consensus_results-notds/ILM3_2024AH06_229E-HTB-46_P3_reads_to_reference_MW202340.1_lofreq.vcf" # submission 955
)

sample_names=(
    "NR-540_vcf_mode" #994
    "HTB-46_vcf_mode" #955
)

references=(
    "AF528167.fasta"
    "MW202340.1.fasta"
)

# for (( i=0; i<${#vcfs[@]}; i++ ));
# do
    
#     vcf="${vcfs[i]}"
#     output_name="${sample_names[i]}"
#     ref="${references[i]}"
#     ref_no_ext="${ref%.*}" 
#     output_dir="$working_dir"/STANDARD/$output_name
#     bash /home/shared/repos/RAISIN_develop/run_raisin.sh \
#         -m "STANDARD" \
#         -s "VCF" \
#         -v $vcf \
#         -o $output_dir \
#         -f "$output_name" \
#         -t 16 \
#         -r /home/shared/references/"$ref_no_ext".fasta \
#         -g /home/shared/references/"$ref_no_ext".gbk

# done 

######################################################?


#########* COMPARE TEST

#! Compare with results in submission 961 and 955 (last sample)
var_txt_1_paths=(
    "/home/shared/service-data/804/roundthree_MW202340.1/ILM3_2023CB06_229E-CCL-185-Passage-1_p3/consensus_results-notds/ILM3_2023CB06_229E-CCL-185-Passage-1_reads_to_reference_MW202340.1_lofreq_variants.txt"
    "/home/shared/service-data/804/roundthree_MW202340.1/ILM3_2023CB07_229E-CRL-3560-Passage-1_p3/consensus_results-notds/ILM3_2023CB07_229E-CRL-3560-Passage-1_reads_to_reference_MW202340.1_lofreq_variants.txt"
    "/home/shared/service-data/843/ILM3_2023CD03_229E-HTB-55_P1_MW202340.1_consensus/consensus_results-notds/ILM3_2023CD03_229E-HTB-55_P1_reads_to_reference_MW202340.1_lofreq_variants.txt"
)


var_txt_2_paths=(
    "/home/shared/service-data/942/ILM3_2024AF05_229E-CCL-185_Passage_3_Day_7_MW202340.1/consensus_results-notds/ILM3_2024AF05_229E-CCL-185_Passage_3_Day_7_reads_to_reference_MW202340.1_lofreq_variants.txt" 
    "/home/shared/service-data/942/ILM3_2024AF06_229E-CRL-3560_Passage_3_Day_7_MW202340.1/consensus_results-notds/ILM3_2024AF06_229E-CRL-3560_Passage_3_Day_7_reads_to_reference_MW202340.1_lofreq_variants.txt"
    "/home/shared/service-data/955/ILM3_2024AH04_229E-HTB-55_P3_MW202340.1/consensus_results-notds/ILM3_2024AH04_229E-HTB-55_P3_reads_to_reference_MW202340.1_lofreq_variants.txt"
)

var_txt_1_names=(
    "Passage_3_CCL-185"
    "Passage_3_CRL-3560"
    "Passage_1_HTB-55"
)

var_txt_2_names=(
    "Passage_3_Day_7_CCL-185"
    "Passage_3_Day_7_CRL-3560"
    "Passage_2_HTB-55"
)


# for (( i=0; i<${#var_txt_1_paths[@]}; i++ ));
# do
    
#     var_txt_1="${var_txt_1_paths[i]}"
#     var_txt_2="${var_txt_2_paths[i]}"
#     var_txt_1_name="${var_txt_1_names[i]}"
#     var_txt_2_name="${var_txt_2_names[i]}"
#     output_name="$var_txt_1_name"_compared_w_"$var_txt_2_name"

#     output_dir="$working_dir"/COMPARE/$output_name
#     ref="${references[i]}"
#     ref_no_ext="${ref%.*}" 
#     bash run_raisin.sh \
#         -m COMPARE \
#         -x $var_txt_1 \
#         -y $var_txt_2 \
#         -o $output_dir \
#         -f $output_name \
#         -X $var_txt_1_name \
#         -Y $var_txt_2_name

# done 

######################################################*

#########! ANCHOR TEST (no download)
ngs_ids=(
    # "ILM3_2024AL06_NR-59711" # 982
    # "ILM3_2024AB07_NR-59705" # 909
    # "ILM3_2023AR07_NR-59317" # 541 # type III variant
    "ILM3_2023AN04_NR-59129" # 524 # type II deletion + insertion
)

references=(
    # "EPI_ISL_18403093.fasta"
    # "EPI_ISL_18432261.fasta"
    # "MT952602.fasta"
    "EPI_ISL_15509864.fasta"
)


for (( i=0; i<${#ngs_ids[@]}; i++ ));
do
    
    ngs_id="${ngs_ids[i]}"
    run_id=$(echo "${ngs_ids[i]}" | cut -c1-11)
    r1_path=/home/shared/raw_fastqs/"$run_id"/"$ngs_id"_R1.fastq.gz
    r2_path=/home/shared/raw_fastqs/"$run_id"/"$ngs_id"_R2.fastq.gz
    anchor_ref="/home/shared/references/NC_045512.2.fasta"
    anchor_gbk="/home/shared/references/NC_045512.2.gbk"
    ref="${references[i]}"
    ref_no_ext="${ref%.*}" 
    output_dir="$working_dir"/ANCHOR/"$ngs_id"
    bash /home/shared/repos/RAISIN_develop/run_raisin.sh \
        -m "ANCHOR" \
        -s "SEQ" \
        -1 $r1_path \
        -2 $r2_path \
        -o $output_dir \
        -f "$ngs_id" \
        -t 16 \
        -r /home/shared/service-data/sars_references/"$ref_no_ext".fasta \
        -a $anchor_ref \
        -k $anchor_gbk


done 

######################################################!