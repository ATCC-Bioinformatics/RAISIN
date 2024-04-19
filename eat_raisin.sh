#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=normal
#SBATCH --job-name=raisin_test

source /opt/conda/etc/profile.d/conda.sh
conda activate oatmeal_2024JAN_V2

ngs_ids=(
    "ILS3_2023AB05_ML-1015" # Flu
)

references=(
    "OP213369.1.fasta"
)


oatmeal_repo=/home/shared/repos/RAISIN_develop
source $oatmeal_repo/oatmeal.sh

cons="/home/shared/service-data/979/ILM3_2024AL03_229E_Calu-3_Passage-2_MW202340.1/consensus_results-notds/ILM3_2024AL03_229E_Calu-3_Passage-2_reads_to_reference_MW202340.1_consensus.fasta"

# cat /home/shared/references/KY684760.fasta > test_mafft_input.fasta
# cat /home/shared/references/MW202340.1.fasta >> test_mafft_input.fasta
# cat $cons >> test_mafft_input.fasta
# mafft test_mafft_input.fasta > test_mafft_alignment.fasta

vcf="/home/shared/service-data/863/ILM3_2023CF07_NR-59685_EPI_ISL_18507690/reads_to_reference/ILM3_2023CF07_NR-59685_reads_to_reference_EPI_ISL_18507690_final.vcf"
mafft="/home/shared/service-data/863/ILM3_2023CF07_NR-59685_EPI_ISL_18507690/reads_to_reference/ILM3_2023CF07_NR-59685_mafft_alignment.fasta"
gbk="/home/shared/references/NC_045512.2.gbk"

outdir="/home/shared/repos/RAISIN_develop"
# python msa_raisin.py
python universal_raisin.py $vcf $gbk "ANCHOR" "full_test" $outdir $mafft


### RAISIN Test Two  #! Try with a different sequence set
# Anchor reference: 
# Strain reference: 

# for (( i=0; i<${#ngs_ids[@]}; i++ ));
# do
    
#     ngs_id="${ngs_ids[i]}"
#     run_id=$(echo "${ngs_ids[i]}" | cut -c1-11)
#     r1_path=/home/shared/raw_fastqs/"$run_id"/"$ngs_id"_R1.fastq.gz
#     r2_path=/home/shared/raw_fastqs/"$run_id"/"$ngs_id"_R2.fastq.gz

#     ref="${references[i]}"
#     ref_no_ext="${ref%.*}" 
#     bash /home/shared/repos/RAISIN_develop/run_raisin.sh \
#         -m "STANDARD" \
#         -s "SEQ" \
#         -1 $r1_path \
#         -2 $r2_path \
#         -o /home/shared/projects/RAISIN/"$ngs_id"_download_ref_with_accs \
#         -f "$ngs_id" \
#         -t 16 \
#         -n "OP213369.1" \
#         -e "dyarmosh@atcc.org" \
#         -d

# done 

# # mafft test_mafft.fasta > test_mafft_output.fasta