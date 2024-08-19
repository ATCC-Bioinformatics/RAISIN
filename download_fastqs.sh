#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=normal
#SBATCH --job-name=download_fq

source /opt/conda/etc/profile.d/conda.sh
conda activate oatmeal_2024JAN_V2

    # "ILM3_2023AR07_NR-59317" # 541 # type III variant
    # "ILM3_2023AN04_NR-59129" # 524 # type II deletion + insertion

python /home/shared/repos/RAISIN_develop/getRawFastqsFromOneCodex.py ILM3_2023AR07_NR-59317 False
python /home/shared/repos/RAISIN_develop/getRawFastqsFromOneCodex.py ILM3_2023AN04_NR-59129 False
