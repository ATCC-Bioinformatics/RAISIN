
.___      .    _   _____ _ __    _
/   \    /|    |  (      | |\   | 
|__-'   /  \   |   `--.  | | \  | 
|  \   /---'\  |      |  | |  \ | 
/   \,'      \ / \___.'  / |   \| 

# Introduction 
ATCC presents RAISIN (Retrieving Amino acid Implications from Sequencing IteratioNs), a simple, fast, and accurate variant annotation pipeline built for characterizing and notating variants given a pair of Illumina sequencing FASTQs and an NCBI reference accession number

# Getting Started
TODO: Guide users through getting your code up and running on their own system. In this section you can talk about:
1.	Installation process
2.	Software dependencies
3.	Latest releases
4.	API references

# Usage
```
Usage: 
    -o for top-level output/working directory (required),
    -f for an output name to add to the beginning of all generated files (optional, default=RAISIN_analysis),
    -t for threads (optional, default=8),

# RAISIN Modes ####################
----------------------------------------------------------------------------------------------
STANDARD mode:
    # Example 1: bash run_raisin.sh -m STANDARD -s SEQ -1 sample1_R1.fastq.gz -2 sample1_R2.fastq.gz -o sample1_results 
                    -f sample1 -n MN02121.1 -e username@gmail.com -d
    # Example 2: bash run_raisin.sh -m STANDARD -s VCF sample1.vcf -o sample1_results 
                    -f sample1 -r MN02121.1.fasta -g MN02121.1.gbk
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
----------------------------------------------------------------------------------------------
ANCHOR mode:
    # Example 1: bash run_raisin.sh -m ANCHOR -s SEQ -1 sample1_R1.fastq.gz -2 sample1_R2.fastq.gz -o sample1_results -f sample1_vs_anchor
                        -n MN02121.1 -b MZ45991.1 -e username@gmail.com -d
    # Example 2: bash run_raisin.sh -m ANCHOR -s SEQ -1 sample1_R1.fastq.gz -2 sample1_R2.fastq.gz -o sample1_results -f sample1_vs_anchor
                        -r MN02121.1.fasta -g MN02121.1.gbk -a MZ45991.1.fasta -k MZ45991.1.gbk
    -m ANCHOR,
    ****** Inputs ******
        -s SEQ (to indicate using sequencing reads),
        -1 for filepath to forward Illumina fastq file,
        -2 for filepath to reverse Illumina fastq file

    ****** For References ******
    Option 1:
        -r for path to strain reference,
        -g for path to strain GBK,
        -a for path to anchor reference,
        -k for path to anchor GBK,
    Option 2:
        -n for NCBI Accession ID to strain reference
        -b for NCBI Accession ID to anchor reference
        -d to download reference and gbk file based on NCBI assession ID,
        -e for Entrez email address to use to download references

----------------------------------------------------------------------------------------------
COMPARE mode:
    # Example 1: bash run_raisin.sh -m COMPARE -x sample1_variants.txt -y sample2_variants.txt 
                -o comp_results -f sample1_vs_sample2 -X Sample_1 -Y Sample_2
    -m standard,
    -x filepath to first variants.txt to compare,
    -y filepath to second variants.txt to compare,
    -X identifier for first variants.txt, i.e. Passage_0
    -Y identifier for second variants.txt, i.e. Passage_1
----------------------------------------------------------------------------------------------

Additional Help:
######## File Inputs  ####################
If user choose either STANDARD or ANCHOR mode, there are two options for file input. Users can choose between either submitting deinterleaved Illumina
sequencing readsets or submitting VCF files. 

######## References ####################
If user choose either STANDARD or ANCHOR mode, a reference and a GenBank file MUST be supplied. A path to both files will suffice BUT
users can also provide an NCBI Accession ID and an Entrez email address instead. RAISIN will then download the reference fasta and
GenBank file based on the given NCBI Accession ID. 

Option 1:
    -r/a for path to reference fasta file (required if not using -d flag),
    -g/k for path to reference gbk file (required if not using -d flag),
Option 2:
    -n/b for NCBI Accession ID, i.e. OP213694,
    -d to download reference and gbk file based on NCBI assession ID or virus, -a or -n flags must be provided,
    -e for Entrez email address to use to download references,

################################
```
# Contribute
TODO: Explain how other users and developers can contribute to make your code better. 

If you want to learn more about creating good readme files then refer the following [guidelines](https://docs.microsoft.com/en-us/azure/devops/repos/git/create-a-readme?view=azure-devops). You can also seek inspiration from the below readme files:
- [ASP.NET Core](https://github.com/aspnet/Home)
- [Visual Studio Code](https://github.com/Microsoft/vscode)
- [Chakra Core](https://github.com/Microsoft/ChakraCore)