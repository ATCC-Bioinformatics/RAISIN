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
    r1_path=/home/shared/raw_fastqs/"$run_id"/"$ngs_id"_R1.fastq.gz
    r2_path=/home/shared/raw_fastqs/"$run_id"/"$ngs_id"_R2.fastq.gz

    ref="${references[i]}"
    ref_no_ext="${ref%.*}" 
    # bash /home/shared/repos/RAISIN_develop/run_raisin.sh \
    #     -1 $r1_path \
    #     -2 $r2_path \
    #     -o /home/shared/projects/RAISIN/"$ngs_id"_test1 \
    #     -f "$ngs_id" \
    #     -r /home/shared/references/"$ref" \
    #     -t 16 \
    #     -g /home/shared/references/$(basename $ref .fasta).gbk
    bash /home/shared/repos/RAISIN_develop/run_raisin.sh \
        -1 $r1_path \
        -2 $r2_path \
        -o /home/shared/projects/RAISIN/"$ngs_id"_download_ref_with_accs \
        -f "$ngs_id" \
        -t 16 \
        -a "OP213369.1" \
        -e "dyarmosh@atcc.org" \
        -d

done 
