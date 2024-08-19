#!/usr/bin/env python
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
# from multiprocessing import Pool
import pandas as pd
import onecodex
import os
import sys
# import pyodbc
from glob import glob
# import sqlalchemy # get Dave to install this

# con = pyodbc.connect('Driver={ODBC Driver 17 for SQL Server};Server=aws-sql1;Database=Pantry;UID=hpcdatareader_svc;PWD=8wZLUxR9O9fo#5RMXJw4;')

ngs_id = sys.argv[1]
download_other_lots = sys.argv[2]   # True/False value
ont_ngs_id = ""

if len(sys.argv) > 3:
    ont_ngs_id = sys.argv[3]

with open('/home/shared/onecodex-api-key/sbc-ocx-api-key.txt') as fr:
    api_key = fr.read()
    
ocx = onecodex.Api(api_key) # Connects to OneCodex API
output_dir = "/home/shared/raw_fastqs"  # all raw data will now go to this folder

incom_mat = pd.read_csv("/home/shared/sbc_scripts/sbc-dashboard/dashboard_cards/Incoming_Materials.txt", sep="\t")
incom_mat = incom_mat[["NUCLEIC_ID", "LOTNUM"]]
df = pd.read_csv("/home/shared/sbc_scripts/genome_status/genome_status.csv")
df = df[["Illumina Sample ID", "ONT Sample ID", "Illumina Sample Filename", "ONT Sample Filename", "External Sample ID"]]

merged_df = df.merge(incom_mat, how="left", left_on="External Sample ID", right_on="NUCLEIC_ID")
merged_df = merged_df[["Illumina Sample ID", "ONT Sample ID", "Illumina Sample Filename", "ONT Sample Filename", "External Sample ID", "LOTNUM"]]

hpc_seq_path="/home/shared/raw_fastqs"
api_samples = pd.read_csv("/home/shared/sbc_scripts/sbc-dashboard/dashboard_cards/onecodex_all_metadata.tsv", sep="\t")

def download_and_deinterleave(lookup_id, ngs_id, output_dir):
    # If we can't find it in genome status, try searching manually through the API by ngs_id

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        print(f"Directory {output_dir} Created.")
    else:    
        print(f"Directory {output_dir} already exists.")

    if lookup_id == "":
        #samples = ocx.Samples.all().filter(lambda x: x.filename.startswith(ngs_id))
        sample = api_samples[api_samples['ngs_id'] == ngs_id]
        sample_lot = sample['lot'].tolist()
        sample_lot = [item for item in sample_lot if item != "nan" ]
        samples = list(set(api_samples[api_samples['lot'].isin(sample_lot)]['ngs_id'].tolist()))
        samples_hex = list(set(api_samples[api_samples['lot'].isin(sample_lot)]['sample_id'].tolist()))
        if len(samples_hex) == 0:
            print(f"Couldn't find {ngs_id} on One Codex. Has this been processed?")
            return
        for i in range(len(samples)):
            if "trimmed" not in samples[i]:
                sample = ocx.Samples.get(samples_hex[i])
                run_id = "_".join(sample.filename.split("_")[0:2])[:-2]
                sample.download(f"{hpc_seq_path}/{run_id}/{sample.filename}")

    else:
        try:
            sample = ocx.Samples.get(lookup_id)
            sample.download(f"{output_dir}/{ngs_id}.fastq.gz")
        except:
            print(f"Ran into issues download via {lookup_id}, trying manually")
            samples = ocx.Samples.all().filter(lambda x: x.filename.startswith(ngs_id))
            if len(samples) == 0:
                print(f"Couldn't find {ngs_id} on One Codex. Has this been processed?")
            for i in range(len(samples)):
                if "trimmed" not in samples[i]:
                    sample = ocx.Samples.get(samples[i].id)
                    sample.download(f"{output_dir}/{ngs_id}.fastq.gz")

    ## Deinterleave fastqs
    if not ngs_id.startswith("ON"):
        for i in range(1,3):
            print("Deinterleaving")
            os.system(f"/home/bin/seqtk seq -{i} {output_dir}/{ngs_id}.fastq.gz > {output_dir}/{ngs_id}_R{i}.fastq")
            os.system(f"chmod 775 {output_dir}/{ngs_id}_R{i}.fastq")
            os.system(f"gzip -f {output_dir}/{ngs_id}_R{i}.fastq")
        
        os.system(f"rm {output_dir}/{ngs_id}.fastq.gz")
    else:
        os.system(f"chmod 775 {output_dir}/{ngs_id}.fastq.gz")
        os.system(f"chmod -R 775 {output_dir}")

# merged_df = pd.read_sql('SELECT gnm_stus.ctlgnum, gnm_stus.nucleic_id, illumina_sample_filnm, illumina_sample_id, ont_sample_filnm, \
#                  ont_sample_id, lotnum FROM GNM_STUS \
#                  LEFT OUTER JOIN INCOMING_MATERIALS ON GNM_STUS.NUCLEIC_ID = INCOMING_MATERIALS.NUCLEIC_ID',con=con)
def download_fastqs(ngs_id, output_dir, download_other_lots):
    run_id = "_".join(ngs_id.split("_")[0:2])[:-2]
    if download_other_lots.lower() == "true":
        if ngs_id.startswith("I"):
            #ngs_info = merged_df[merged_df["Illumina Sample Filename"] == f"{ngs_id}.fastq.gz"]
            ngs_info = api_samples[api_samples['ngs_id'] == ngs_id]
        else:
            #ngs_info = merged_df[merged_df["ONT Sample Filename"] == f"{ngs_id}.fastq.gz"]
            ngs_info = api_samples[api_samples['ngs_id'] == ngs_id]
        if len(ngs_info) > 0:
            lotnum = ngs_info["lot"].head(1).item()  # Grab the first entry
            if str(lotnum) != "nan":
                matching_lots = api_samples[api_samples["lot"] == lotnum]
                print(f"Here are the samples associated with {ngs_id}; LotNo. {lotnum} on OneCodex")
                print(matching_lots)
                print()
                # Retrieve all fastqs for that lotnum and deduplicate list
                #! Added to fix doing this twice for ONT and ILM
                if ngs_id.startswith("I"):
                    fastqs = list(set(matching_lots[matching_lots['platform'] == "Illumina"]['ngs_id'].tolist()))
                else:
                    fastqs = list(set(matching_lots[matching_lots['platform'] == "ONT"]['ngs_id'].tolist()))
            else:
                fastqs = [f"{ngs_id}.fastq.gz"]
            
            for file in fastqs:
                if str(file) == "nan":
                    continue
                print("looping through ngsid: ", file)
                run_id = "_".join(file.split("_")[0:2])[:-2]
                if len(run_id) < 10:
                    run_id = "_".join(file.split("_")[0:2])
                current_ngs_id = file.split(".")[0]
                # Check if file is on the HPC first
                hpc_file_path = glob(f"{hpc_seq_path}/{run_id}/{current_ngs_id}*.fastq.gz")
                #raw_path = glob(f"{output_dir}/{current_ngs_id}*.fastq.gz")
                if len(hpc_file_path) == 0: # and len(raw_path) == 0: #! if not on hpc or in /home/shared/raw_fastqs
                    # If not on hpc, check given output dir next
                    print(f"Files not found on HPC, searching in {output_dir}")
                    #files_output_path = glob(f"{output_dir}/{current_ngs_id}*.fastq.gz")
                    # if it still can't be found, download from OC
                    #if len(files_output_path) == 0:
                    if str(lotnum) != "nan":
                        if current_ngs_id.startswith("I"):
                            sample_id = matching_lots[matching_lots['ngs_id'] == current_ngs_id]['sample_id'].iloc[0]
                        else:
                            sample_id = matching_lots[matching_lots['ngs_id'] == current_ngs_id]['sample_id'].iloc[0]
                    else:
                        if current_ngs_id.startswith("I"):
                            sample_id = api_samples[api_samples["ngs_id"] == current_ngs_id]['sample_id'].iloc[0]
                            #sample_id = matching_lots[matching_lots['ngs_id'] == current_ngs_id]['sample_id'].iloc[0]
                        else:
                            sample_id = api_samples[api_samples["ngs_id"] == current_ngs_id]['sample_id'].iloc[0]

                    print(f"Downloading fastqs for {current_ngs_id} to {output_dir}/{run_id}")
                    download_and_deinterleave(sample_id, current_ngs_id, f"{output_dir}/{run_id}")
                    # elif len(raw_path) > 0:
                    #     print(f"Fastqs for {current_ngs_id} already in {output_dir}")
                    # else:
                    #     print(f"Fastqs for {current_ngs_id} already in {hpc_seq_path}")

                else:
                    print(f"Fastqs for {current_ngs_id} already on HPC")
        else:
            print(f"Whoops, couldn't find {ngs_id} in the genome status. Searching for it on One Codex the long way")
            sample_id = ""
            download_and_deinterleave(sample_id, ngs_id, f"{output_dir}/{run_id}")

    else:     
        # Since searching through the API takes way too long, it's much quicker to get the sample ID from gnm_stus
        if ngs_id.startswith("IL"):
            sample_df = df[df["Illumina Sample Filename"].str.contains(ngs_id, na=False)]  # ignore NA values with na=False
            lookup_id = "".join(sample_df["Illumina Sample ID"].drop_duplicates())
        elif ngs_id.startswith("ON"):
            sample_df = df[df["ONT Sample Filename"].str.contains(ngs_id, na=False)]  # ignore NA values with na=False
            lookup_id = "".join(sample_df["ONT Sample ID"].drop_duplicates())
        
        download_and_deinterleave(lookup_id, ngs_id, f"{output_dir}/{run_id}")


# Download Illumina fastqs
download_fastqs(ngs_id, output_dir, download_other_lots)

if len(ont_ngs_id) > 3: # download ONT only if ont ngs id is given
    download_fastqs(ont_ngs_id, output_dir, download_other_lots)


