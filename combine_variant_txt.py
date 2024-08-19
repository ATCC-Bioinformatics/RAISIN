import pandas as pd
from glob import glob
from functools import reduce
import os
import sys

# Read in arguments
variant_txt_1 = sys.argv[1]
variant_txt_2 = sys.argv[2]
variant_id_1 = sys.argv[3]
variant_id_2 = sys.argv[4]
outdir = sys.argv[5]

variant_dataframes = []
variant_txts = {variant_txt_1: variant_id_1, variant_txt_2: variant_id_2}
for txt, var_id in variant_txts.items():
    variant_df = pd.read_csv(txt, sep="\t")
    variant_df = variant_df.rename(columns={"Allele frequency": f"{var_id}_Allele frequency", "AA mutation": f"{var_id}_AA Mutation"})
    variant_dataframes.append(variant_df)

df_merged = reduce(lambda left,right: pd.merge(left,right,on=['Position','Reference allele','Alternate allele','Protein'],
            how='outer'), variant_dataframes)

# Reorder columns
column_to_move = df_merged.pop("Protein")
df_merged.insert(3, "Protein", column_to_move)
# Calculate frequency difference
df_merged[f"{variant_id_1}-{variant_id_2} Frequency Difference"] = df_merged[f"{variant_id_1}_Allele frequency"] - df_merged[f"{variant_id_2}_Allele frequency"]

df_merged.to_csv(f"{outdir}/{variant_id_1}-{variant_id_2}_combined_variant_tables.tsv", index=False, sep="\t")
