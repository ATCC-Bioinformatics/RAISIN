# The RAISIN Pipeline
Retrieving Amino acid Implications from Sequencing IteratioNs (RAISIN): A Pipeline Intended to Better Characterize Viral Genomic Variants

ATCC presents RAISIN (Retrieving Amino acid Implications from Sequencing IteratioNs), a simple, fast, and accurate variant annotation pipeline built for characterizing and notating viral genomic variants. As it stands, RAISIN is currently meant for all **viruses**, but we plan to add support for bacteria and other microbes in the future. 

```
.___      .    _   _____ _ __    _
/   \    /|    |  (      | |\   | 
|__-'   /  \   |   `--.  | | \  | 
|  \   /---'\  |      |  | |  \ | 
/   \,'      \ / \___.'  / |   \| 
```
# Overview

RAISIN has three modes: STANDARD, ANCHOR, and COMPARE, each of which allow users to get a different glimpse of the variants for a sample.

![RAISIN Overview](https://github.com/ATCC-Bioinformatics/RAISIN/blob/develop/readme_images/RAISIN%20OVERVIEW_v3.jpg)

# Installation
RAISIN is pipeline with mainly Bash and Python scripts. We recommend installing all dependencies through conda or mamba (for faster installation).

First, use conda to create your conda environment:
```
conda create --name raisin_env
conda activate raisin_env
```
Then, install the following packages in the environment. We recommend installing with mamba for faster installation. 
```
mamba install gatk4
mamba install -c bioconda lofreq bwa bedtools qualimap fastp bbmap picard samtools mafft fastqc bcftools
mamba install multiqc
mamba install -c conda-forge biopython pandas
mamba install conda-forge::r-gplots
mamba install conda-forge::r-gsalib
```
# Usage
### STANDARD mode
In STANDARD mode, users can either use Illumina sequencing FASTQs or a VCF as the main input. If given the reads, STANDARD RAISIN will map the reads to a user-given reference, call variants, and then call consensus. The generated VCF is then analyzed against the user-given annotation file to better characterize each variant. If a VCF is given, the consensus step is skipped, and the mode moves straight to variant analysis. The main output of STANDARD mode is a TSV containing all variants with their respective genomic region and amino acid translation. 

![RAISIN STANDARD Mode](https://github.com/ATCC-Bioinformatics/RAISIN/blob/develop/readme_images/STANDARD_V5.jpg)

STANDARD mode is great if you're looking to obtain additional metadata for your viral variants, specifically metadata on where variants occur in the genome and the resulting amino acid mutations. 

#### Examples: 
Example if choosing to run with Illumina readset and would like to download reference + annotation from NCBI: 
```
bash raisin.sh -m STANDARD -s SEQ -1 sample1_R1.fastq.gz -2 sample1_R2.fastq.gz -o sample1_results \
                    -f sample1 -n MN02121.1 -e username@gmail.com -d
```
Example if choosing to run with vcf: 
```
bash raisin.sh -m STANDARD -s VCF -v sample1.vcf -o sample1_results \
                    -f sample1 -r MN02121.1.fasta -g MN02121.1.gbk
```

#### Available arguments:
```
    -o for top-level output/working directory (required),
    -f for an output name to add to the beginning of all generated files (optional, default=RAISIN_analysis),
    -t for threads (optional, default=8),
    -q for minimum variant frequency threshold, enter as a decimal, i.e. 0.05, (optional, default=0.05),
    -c for minimum read coverage for variant calling, (optional, default=10) 
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
```

### COMPARE mode
In COMPARE mode, variants from two samples are compared together. This mode requires two RAISIN STANDARD mode output TSVs as input so we recommend running STANDARD mode for both samples first and then running COMPARE afterwards. Any samples that are compared with COMPARE mode must use the same reference when variants were called in STANDARD mode. In other words, if you were to compare sample 1 and sample 2, they must both have variants called using the same reference. The main output of COMPARE mode is a TSV with an additional column displaying the frequncy difference between the two samples.

![RAISIN COMPARE Mode](https://github.com/ATCC-Bioinformatics/RAISIN/blob/develop/readme_images/COMPARE_V4.jpg)

COMPARE mode is ideal in case where you are comparing two passages or samples against each other to not only determine what variants are shared with each other, but to see the frequency differences of those variants between two samples and to see if any variants were lost or gained between the two samples. 

#### Example: 
```
bash raisin.sh -m COMPARE -x sample1_variants.txt -y sample2_variants.txt \
                    -o comp_results -f sample1_vs_sample2 -X Sample_1 -Y Sample_2
```
#### Available arguments:

```
    -o for top-level output/working directory (required),
    -f for an output name to add to the beginning of all generated files (optional, default=RAISIN_analysis),
    -t for threads (optional, default=8),
    -m COMPARE,
    -x filepath to first variants.txt to compare,
    -y filepath to second variants.txt to compare,
    -X identifier for first variants.txt, i.e. Passage_0
    -Y identifier for second variants.txt, i.e. Passage_1

```
### ANCHOR mode
In ANCHOR mode, a three-way comparison of the anchor sequence, strain sequence, and consensus sequence is performed. Inputs for anchor mode are the path to Illumina sequencing readset, the path to the strain reference, the path to anchor reference, and the path to anchor GBK. 

![RAISIN ANCHOR Mode](https://github.com/ATCC-Bioinformatics/RAISIN/blob/develop/readme_images/ANCHOR_V4.jpg)

ANCHOR mode is a bit of special case but it is ideal to use in situations where you want to compare your sample to a strain reference and a species reference. A strain reference would be a sequence that is the most similar to your sample. A species reference would be a sequence that is most similar to the expected species of your sample. With these sequences, ANCHOR mode will return the variants from each. The resulting output file will show the variants found and will report all variants using the coordinate sequence of the anchor/species reference.

ANCHOR mode uses a variant classification system (I-IV) to differentiate between variants. For definition of each classification, see below: 

![RAISIN ANCHOR Mode](https://github.com/ATCC-Bioinformatics/RAISIN/blob/develop/readme_images/variant_classifications_anchor.png)

*Legend:*
* **ANCHOR**: a scientific community recognized or experiment-specific sequence
* **STRAIN**: a sequence that is more closely related to the sample
* **SAMPLE**: a consensus sequence of the sample reads mapped to the strain sequence

#### Examples: 
Example to run ANCHOR mode and download reference + annotation from NCBI: 
```
bash raisin.sh -m ANCHOR -s SEQ -1 sample1_R1.fastq.gz -2 sample1_R2.fastq.gz -o sample1_results -f sample1_vs_anchor \
                          -n MN02121.1 -b MZ45991.1 -e username@gmail.com -d
```
Example to run ANCHOR mode and use filepaths for reference + annotation:
```
bash raisin.sh -m ANCHOR -s SEQ -1 sample1_R1.fastq.gz -2 sample1_R2.fastq.gz -o sample1_results -f sample1_vs_anchor \
                          -r MN02121.1.fasta -a MZ45991.1.fasta -k MZ45991.1.gbk
```

#### Available arguments:
```
    -o for top-level output/working directory (required),
    -f for an output name to add to the beginning of all generated files (optional, default=RAISIN_analysis),
    -t for threads (optional, default=8),
    -q for minimum variant frequency threshold, enter as a decimal, i.e. 0.05, (optional, default=0.05),
    -c for minimum read coverage for variant calling, (optional, default=10) 
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

  ```

# Additional Help:

#### References and Annotations
    
If user choose either STANDARD or ANCHOR mode, a reference and a GenBank file MUST be supplied. A path to both files will suffice BUT users can also provide an NCBI Accession ID and an Entrez email address instead. RAISIN will then download the reference fasta and GenBank file based on the given NCBI Accession ID. 
  
Option 1:
  -r/a for path to reference fasta file (required if not using -d flag), (use -a for reference path for anchor)
  -g/k for path to reference gbk file (required if not using -d flag), (use -k for annotations path for anchor)
Option 2:
  -n/b for NCBI Accession ID, i.e. OP213694,
  -d to download reference and gbk file based on NCBI assession ID or virus, -a or -n flags must be provided,
  -e for Entrez email address to use to download references,
