from Bio import Entrez
import pandas as pd
import sys
import os

acc = sys.argv[1]
email = sys.argv[2]
outdir = sys.argv[3]
log_path = sys.argv[4]
Entrez.email = email

if not os.path.exists(f'{outdir}/{acc}.fasta') and not os.path.exists(f'{outdir}/{acc}.gbk'):
  try:
      id = Entrez.read(Entrez.esearch(db='nucleotide',term=acc))['IdList'][0]
      os.system(f"echo 'Retrieving fasta and GBK files...' >> {log_path}")
      fasta = Entrez.efetch(db='nucleotide',id=id,retmode='text',rettype='fasta').read()
      gb = Entrez.efetch(db='nucleotide',id=id,retmode='text',rettype='gbk').read()
      with open(f'{outdir}/{acc}.fasta','w') as out:
        out.write(fasta)
      with open(f'{outdir}/{acc}.gbk','w') as out:
        out.write(gb)
  except IndexError:
      os.system(f"echo 'Provided accession number is invalid: {acc}' >> {log_path}")
