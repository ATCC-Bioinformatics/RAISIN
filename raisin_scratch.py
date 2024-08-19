# Python script to generate variant_table.csv for BEI reports.
# It takes in the path to the MSA output and determines variant info
# by iterating down the alignment and comparing the nucleotides at each
# position in each sequence. The nt info is then compared to data from
# Marco's table (saved in BEI_sars_genes.pkl) to determine AA changes.

# Aligned sequences are truncated from 1-50, and 29800-end due to strain
# and consensus typically being shorter than anchor sequence.
import pickle as pkl
import argparse
from Bio import pairwise2
# from Bio.SubsMat import MatrixInfo as matlist
# matrix = matlist.blosum62

# SARS-COV-2 data from Marco's table ------all values are strings
# sars_genes[position] = {'Gene':___,
#                    'Region':___,
#                    'nt':___,
#                    'Codon Position':___,
#                    'Codon':___,
#                    'AA':___,
#                    'AA #':___}

BEI_file = open('/home/shared/BEI_Pipeline/BEI_sars_genes_2.pkl', 'rb')
# positions start at 1
sars_genes = pkl.load(BEI_file)

## Genetic Code ----------------------------------------
genetic_code = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }

# iupac_code = {
#   'M':'AC',
#   'R':'AG',
#   'W': 'AT',
#   'S':'CG',
#   'Y':'CT',
#   'K':'GT',
#   'V':'ACG',
#   'H':'ACT',
#   'D':'AGT',
#   'B':'CGT',
#   'N':'ACGT'
# }

def get_variants(vcf_path):
    ### Collect freebayes variants, allele freq, and depth
    try:
        # read in file
        lofreq = open(vcf_path,"r")
        variants = {}
        # bypass comments at beginning of file
        for line in lofreq:
            if '#' == line[0]:
                pass
            else:
                # create variables based upon column information
                line = line.split('\t')
                pos = int(line[1])
                ref = line[3].lower()
                alt = line[4].lower()
                cov = line[7].split(';')[0][3:] # 2: to remove DP=
                af = line[7].split(';')[1][3:] # 2: to remvove AF=
                # change significant figures to 4
                if af == '1.0':
                    af = '1.0000'
                if pos in variants.keys():
                    # convert to list of lists to contain all variants at that position
                    if isinstance(variants[pos][0],list):
                        variants[pos].append([pos,cov,af,ref,alt])
                    else:
                        variants[pos] = [variants[pos]]
                        variants[pos].append([pos,cov,af,ref,alt])
                else:
                    variants[pos] = [pos,cov,af,ref,alt]
        # print(variants)
        lofreq.close()
        return variants
    except:
        print("Please provide a valid lofreq file path.")
        return ''



# Take in the MSA and create a dict of the three alignment seqs (anchor, strain, and consensus)
def process_alignment(msa_path):
    alignments = {}
    c = 0
    with open(msa_path,'r') as f:
        for line in f:
            # transform nt sequences to one-line, i.e. remove header and carriage returns
            if '>' in line:
                alignments[c] = ''
                c+=1
            else:
                alignments[c-1] += line.lower().replace('\n','')
    return alignments[0],alignments[1],alignments[2]

# Use the sars_genes object to determine gene region, codon, and
# mutation information
# def determine_codon(position,sample_nt):
#     # gather and format gene region info
#     if 'Intergenic' in sars_genes[position]['Gene']:
#         gene_region = sars_genes[position]['Gene'] + sars_genes[position]['Region']
#     elif sars_genes[position]['Region'] != '':
#         gene_region = sars_genes[position]['Gene'] + '(' \
#                     + sars_genes[position]['Region'] + ')'
#     else:
#         gene_region = sars_genes[position]['Gene']
#     # deletion
#     if sample_nt == '-':
#         return gene_region, 'n/a'
#     else:
#         # grab Codon
#         codon = sars_genes[position]['Codon']
#         # check if position is in gene
#         if codon != '':
#             mutated_codon = ''
#             # loop through codon position (1, 2 , 3)
#             for i in range(len(codon)):
#                 # if at the mutated position, then add the sample nt
#                 if i == int(sars_genes[position]['Codon Position']) -1:
#                     mutated_codon += sample_nt
#                 # else, add the sars nt
#                 else:
#                     mutated_codon += codon[i]
#             # translate anchor codon
#             AA = genetic_code[codon.upper()]
#             if 'N' in mutated_codon.upper():
#                 return gene_region, 'ambig'
#             else:
#                 # translate mutated/consensus codon
#                 mutated_AA =  genetic_code[mutated_codon.upper()]
#                 if AA == mutated_AA:
#                     return gene_region, '--'
#                 else:
#                     AA_mutation = AA + sars_genes[position]['AA #'] + mutated_AA # Example output notation: A12Y
#                     return gene_region, AA_mutation
#         else:
#             return gene_region, 'n/a'

# def determine_ins_aa_mutation(insertion_nucs,insertion_position):
#     if 'N' in insertion_nucs:
#         return 'Ambiguous bases'
#     print("insertion position", insertion_position)
#     print("Data:", sars_genes[str(insertion_position)])
#     # If intergenic, then no aa mutation
#     if 'Intergenic' in sars_genes[str(insertion_position)]['Gene'] or 'UTR' in sars_genes[str(insertion_position)]['Gene']:
#         return 'n/a'
#     # if the insertion is a multiple of 3 and occurs immediately after a codon,
#     # then the notation is just a translation of the insertion_nucs
#     if len(insertion_nucs) % 3 == 0 and sars_genes[str(insertion_position)]['Codon Position'] == '3':
#         codons = [''.join(insertion_nucs[i:i+3]) for i in range(0,len(insertion_nucs),3)]
#         aas = ''.join([genetic_code[e] for e in codons])
#         if '*' not in aas:
#             final_notation = sars_genes[str(insertion_position)]['AA #'] + '[' + aas +']'+ sars_genes[str(insertion_position+1)]['AA #']
#         else:
#             final_notation = 'Manually inspect'
#     # if the insertion is a multiple of 3 and occurs inside of a codon, then that codon is mutated
#     # and a string of new codons is inserted. Notation e.g.,  L11V 11[KVR]12
#     elif len(insertion_nucs) % 3 == 0 and sars_genes[str(insertion_position)]['Codon Position'] != '3':
#         # Generate the mutated nucleotide sequence. If the insertion occurs in a codon, then the mutated
#         # nucleotide sequence = original codon up to the insertion + the insertion + the original codon after the insertion

#         # INS between codon position 2 and 3
#         if sars_genes[str(insertion_position)]['Codon Position'] == '2':
#             previous_1_codon_position = insertion_position - 1
#             mutated_nt_sequence = sars_genes[str(insertion_position-1)]['nt'] + sars_genes[str(insertion_position)]['nt'] + insertion_nucs + sars_genes[str(insertion_position+1)]['nt']
#         # INS between codon position 1 and 2
#         else:
#             previous_1_codon_position = insertion_position
#             mutated_nt_sequence = sars_genes[str(insertion_position)]['nt'] + insertion_nucs + sars_genes[str(insertion_position+1)]['nt'] + sars_genes[str(insertion_position+2)]['nt']
#         # Convert mutated nucleotide sequence to amino acid sequence
#         mutated_nt_sequence = mutated_nt_sequence.upper()
#         print("mutated nt sequence: ", mutated_nt_sequence)
#         mutated_codons = [''.join(mutated_nt_sequence[i:i+3]) for i in range(0,len(mutated_nt_sequence),3)]
#         mutated_aa = ''.join([genetic_code[e] for e in mutated_codons])
#         print("mutated_aa sequence: ", mutated_aa)
#         # Generate notation
#         # LV 11-12 KP  11[LVL]12
#         original_pre_aa = sars_genes[str(previous_1_codon_position)]['AA']
#         original_pre_aa_num = sars_genes[str(previous_1_codon_position)]['AA #']
#         final_notation = ''
#         if str(original_pre_aa) != mutated_aa[0]:
#             final_notation = str(original_pre_aa) + str(original_pre_aa_num) + mutated_aa[0] + ' '
#         final_notation += original_pre_aa_num + '[' + mutated_aa[1:] + ']' + str(int(original_pre_aa_num)+1)
#     # if the insertion is not a multiple of 3
#     else:
#         # Find the starting position of the gene in which the insertion occurs
#         position = insertion_position
#         not_found = True
#         while not_found == True:
#             codon_position = sars_genes[str(position)]['Codon Position']
#             aa_num = sars_genes[str(position)]['AA #']
#             if aa_num == '1' and codon_position == '1':
#                 gene_start_position = position
#                 not_found = False
#             else:
#                 # if aa_num != 1 and codon_position != go back a
#                 # position until aa_num == 1 and codon_position == 1
#                 position -= 1
#         # Determine the original nucleotide and amino acid sequence for that gene
#         original_nt_sequence = ''
#         original_aa_sequence = ''
#         stop = False
#         position = gene_start_position
#         print('building original nt and aa seqs')
#         while stop == False:
#             original_nt_sequence += sars_genes[str(position)]['nt']
#             if sars_genes[str(position)]['Codon Position'] == '1':
#                 # Not all genes end with a stop codon. Look for AA # 1 in these instances to determine the end of the gene
#                 if position > gene_start_position and sars_genes[str(position)]['AA #'] == '1':
#                     stop = True
#                 else:
#                     original_aa_sequence += sars_genes[str(position)]['AA']
#             position +=1
#         # Determine the mutated nucleotide sequence
#         mutated_nt_sequence = ''
#         mutated_aa_sequence = ''
#         stop = False
#         position = gene_start_position
#         while stop == False:
#             # once we have gotten to the insertion position, we need to add the
#             # insertion nucleotides AND the nucleotide at insertion_position + 1
#             if position == insertion_position + 1:
#                 # adds insertion nucleotides and the nucleotide directly after the insertion from the original sequence
#                 mutated_nt_sequence += insertion_nucs + sars_genes[str(position)]['nt']
#             elif position < 29800:
#                 # add all nucleotides before and after insertion
#                 mutated_nt_sequence += sars_genes[str(position)]['nt']
#             else:
#                 stop = True # stop at end of anchor sequence
#             position += 1
#         # Convert mutated nucleotide sequence to amino acid sequence
#         for i in range(0,len(mutated_nt_sequence),3):
#             aa = genetic_code[mutated_nt_sequence[i:i+3].upper()]
#             mutated_aa_sequence += aa
#             # stop translating once a stop codon is reached in mutated_aa_sequence
#             if aa == '*':
#                 mutated_nt_sequence = mutated_nt_sequence[:i+3]
#                 break

#         # No notation determined, so just output the entire original and mutated aa sequences
#         final_notation = 'Manually determine. Original aa sequence: '+ original_aa_sequence #### ask Marco
#         final_notation += ' Mutated aa sequence: '+ mutated_aa_sequence #### ask Marco
#     print("final notation: ",final_notation)
#     return final_notation

# def determine_del_aa_mutation(start_position,length):
#     # If intergenic, then no aa mutation
#     if 'Intergenic' in sars_genes[str(start_position)]['Gene'] or 'UTR' in sars_genes[str(start_position)]['Gene']:
#         return 'n/a'
#     not_found = True
#     position = start_position
#     # find gene start position of gene where deletion occurs
#     while not_found == True:
#         codon_position = sars_genes[str(position)]['Codon Position']
#         aa_num = sars_genes[str(position)]['AA #']
#         if aa_num == '1' and codon_position == '1':
#             gene_start_position = position
#             not_found = False
#         else:
#             position -= 1
#     original_nt_sequence = ''
#     original_aa_sequence = ''
#     stop = False
#     position = gene_start_position
#     # determine original nucleotide sequence
#     while stop == False:
#         original_nt_sequence += sars_genes[str(position)]['nt']
#         if sars_genes[str(position)]['Codon Position'] == '1':
#             if position>gene_start_position and sars_genes[str(position)]['AA #'] == '1':
#                 stop = True
#             else:
#                 original_aa_sequence += sars_genes[str(position)]['AA']
#         position +=1

#     mutated_nt_sequence = ''
#     mutated_aa_sequence = ''
#     stop = False
#     position = gene_start_position
#     # determine mutated nucleotide sequence and translated aa sequence
#     while stop == False:
#         if position < start_position or position > start_position+length-1:
#             mutated_nt_sequence += sars_genes[str(position)]['nt']
#             if len(mutated_nt_sequence)%3 == 0:
#                 # not all gene regions end in a stop codon so we check for the next AA # of 1 as a proxy for stop
#                 if position > gene_start_position+3 and sars_genes[str(position)]['AA #'] == '1':
#                     stop = True
#                     # stop creating mutated_aa_sequence when nts translate to stop codon
#                 elif genetic_code[mutated_nt_sequence[-3:].upper()] == '*':
#                     mutated_aa_sequence += genetic_code[mutated_nt_sequence[-3:].upper()]
#                     stop = True
#                 else:
#                     # translate last 3 nts in sequence to get aa, then add to mutated_aa_sequence
#                     mutated_aa_sequence += genetic_code[mutated_nt_sequence[-3:].upper()]
#         elif position > 29799:
#             stop = True
#         position += 1

#     # if the length of the deletion is a multiple of 3, align the original_aa_sequence and mutated_aa_sequence up until the stop codon of each
#     if int(length)%3 == 0:
#         # -.5 gap penalty
#         # -.1 extend penalty
#         # alignment = pairwise2.align.globalds(''.join([e for e in original_aa_sequence if e != '*']),''.join([e for e in mutated_aa_sequence if e != '*']), matrix,-1, -.1,one_alignment_only=True)[0]
#         alignment = pairwise2.align.globalms(''.join([e for e in original_aa_sequence if e != '*']),''.join([e for e in mutated_aa_sequence if e != '*']), 2,-.5,-1, -.1,one_alignment_only=True)[0]
#         original_aligned = alignment[0]
#         mutated_aligned = alignment[1]
#         # print(original_aligned)
#         # print(mutated_aligned)
#         original2mutated = []
#         deleted = []
# #        insertion = []
#         for i in range(len(original_aligned)):
#             # if the nucleotide positions are equal, keep going
#             if original_aligned[i] == mutated_aligned[i]:
#                 pass
# #            elif original_aligned[i] == '-':   #### do we need this code here? ####
# #                insertion.append([mutated_aligned[i],str(i+1)])
#             elif mutated_aligned[i] == '-':
#                 # add to the deleted list a list of aa at i position, and the positional information
#                 deleted.append([original_aligned[i],str(i+1)])
#             else:
#                 # if original aa and mutated aa are not the same append aa of original at i position,
#                 # positional information, and aa of mutated at i position
#                 original2mutated.append([original_aligned[i],str(i+1),mutated_aligned[i]])
#         deletion_notation = ''
#         if len(deleted) > 0:
#             if len(deleted) == 1:
#                 deletion_notation += '∆' + deleted[0][0] + deleted[0][1] # ∆Y144 - deleted aa and deleted aa # if only 1 aa deletion
#             else:
#                 deletion_notation += '∆' + ''.join([e[0] for e in deleted])
#                 deletion_notation += ' (aa' + deleted[0][1] + '-' + deleted[-1][1] + ')' # ∆ABCD (aa1-4) - aa deleted, and aa positions
#         substitution_notation = ''
#         if len(original2mutated) > 0:
#             if len(original2mutated) == 1:
#                 substitution_notation += ''.join(original2mutated[0]) # Y144H
#             else:
#                 substitution_notation += ''.join([e[0] for e in original2mutated])
#                 substitution_notation += ' '+ original2mutated[0][1] + '-' + original2mutated[-1][1] + ' '
#                 substitution_notation += ''.join([e[2] for e in original2mutated]) # ABCDEFG 1-7 HIJKLMN
#         final_notation = substitution_notation + ' ' + deletion_notation # ABCDEFG 1-7 HIJKLMN   ∆ABCD (aa1-4)
#     else:
#         # if length of deletion not divisible by 3
#         different = False
#         i = 0
#         # determine first aa position where original and mutated are not the same
#         while(different==False):
#             if i > len(original_aa_sequence) or i > len(mutated_aa_sequence):
#                 different=True
#             elif original_aa_sequence[i] == mutated_aa_sequence[i]:
#                 i+=1
#                 if i == len(original_aa_sequence):
#                     final_notation = 'Manually determine'
#                     different=True
#             else:
#                 different=True
#         first_different_aa_position = i
#         # original_aa_sequence=original_aa_sequence[i:]
#         # mutated_aa_sequence=mutated_aa_sequence[i:]
#         sequence_notation = ''
#         del_notation = ''
#         # Situation where the mutation causes an early stop codon
#         # mutated sequence in shorter than original sequence
#         if len(original_aa_sequence[first_different_aa_position:]) > len(mutated_aa_sequence[first_different_aa_position:]):
#             mutated_sequence_length = len(mutated_aa_sequence[first_different_aa_position:])-1
#             sequence_notation = str(first_different_aa_position+1) # position of first different aa, e.g. 10
#             sequence_notation += '-' + str(first_different_aa_position+mutated_sequence_length+1) # 10 - 15 adds positions through length of mutated sequence length
#             del_notation = '∆aa' + str(first_different_aa_position+mutated_sequence_length+2)
#             del_notation += '-' + str(len(original_aa_sequence)) # ∆aa16-20 aa positions deleted from original including stop codon through the end of the gene
#             # ABCDEFG original amino acid sequence up to the aa position in the mutated sequence that creates a stop codon
#             final_notation = original_aa_sequence[first_different_aa_position:first_different_aa_position+mutated_sequence_length+1] + ' '
#             final_notation += sequence_notation + ' ' + mutated_aa_sequence[first_different_aa_position:] + ' ' + del_notation # ABCDE 10-15 HIJK* ∆aa16-20
#             # if mutated aa sequence is longer than the original aa sequence, determine notation manually
#         elif  len(mutated_aa_sequence[first_different_aa_position:]) > len(original_aa_sequence[first_different_aa_position:]):
#             final_notation = 'Manually determine'  ##### ask Marco
#     return final_notation

# generate the data for the table
def generate_table(anchor,strain,consensus,variants,low_cov_n_regions):
    out_lines = []  # list of lines that will be output
    apos = 0        # anchor position
    spos = 0        # strain position
    cpos = 0        # consensus position

    for i in range(len(anchor)):
        # NTs
        a = anchor[i]
        s = strain[i]
        c = consensus[i]
        # Incrememnt position
        if a != '-':
            apos+=1
        if s != '-':
            spos+=1
        if c != '-':
            cpos+=1

        coverage = 'n/a'
        sample_frequency = 'n/a'

        # Determine sample frequency and coverage
        for key in variants.keys():
            # variants from lofreq final.vcf, variants[pos] = [pos,cov,af,ref,alt]
            variant = variants[key]
            # check for multiple variants at position
            if isinstance(variant[0],list):
                vpos = variant[0][0]
            else:
                vpos = variant[0]

            ref_allele = variant[-2]
            alt_allele = variant[-1]
            
        # if cons = ref = anchor
        if c == r and r == w:
            # Not a variant
            continue
        if c == 'n' and spos not in list(variants.keys()):
            continue
        if apos == 0:
            continue
        # if ref = n and cons does not = r
        elif r == 'n' and c != 'n':
            # Type N, modified
            if w == '-':
                # determine the gene region of the insertion
                gene_region, AA_mutation = determine_codon(str(apos),c)
                out_line = '\t'.join(['N','INS',str(apos)+c,coverage,'1',sample_frequency,gene_region,'TBD',modified])
                out_lines.append(out_line)
            elif c == '-':
                # determine the gene region of the deletion
                gene_region, AA_mutation = determine_codon(str(apos),c)
                out_line = '\t'.join(['N','DEL',str(apos),coverage,'1',sample_frequency,gene_region,'TBD',modified])
                out_lines.append(out_line)
            else:
                # determine SNP gene region and aa mutation, if w != c and r = n
                variant = w + str(apos) + c
                gene_region, AA_mutation = determine_codon(str(apos),c)
                out_line = '\t'.join(['N','SNP',variant,coverage,'1',sample_frequency,gene_region,AA_mutation,modified])
                out_lines.append(out_line)
        # if ref differs from cons and anchor
        elif c != r and c == w :
            # Type III, modified
            if c == '-':
              # Deal with deletions and insertions later
              out_line = '\t'.join(['III','DEL','TBD','TBD','TBD','TBD','TBD','TBD',modified])
              out_lines.append(out_line)
            #   varNum='III'
            elif r == '-':
              # Deal with deletions and insertions later
              out_line = '\t'.join(['III','INS','TBD','TBD','TBD','TBD','TBD','TBD',modified])
              out_lines.append(out_line)
            #   varNum='III'
            else:
                # if variant is SNP or SUB
                sample_frequency = variants[spos][2]
                coverage = variants[spos][1]
                variant = w + str(apos) + r + '_rev_' + c + '‡' # a1284t_rev_a‡
                gene_region, AA_mutation = determine_codon(str(apos),c)
                out_line = '\t'.join(['III','SNP',variant,coverage,'1',sample_frequency,gene_region,AA_mutation,modified])
                out_lines.append(out_line)
        # if none of the sequences have the same nt at a given position
        elif c != r and c != w and r != w:
            # Type IV, modified
            if c == '-':
              # Deal with deletions and isertions later
              out_line = '\t'.join(['IV','DEL',str(apos),'TBD','TBD','TBD','TBD','TBD',modified])
              out_lines.append(out_line)
            #   varNum='IV'
            elif r == '-':
              # Deal with deletions and isertions later
              out_line = '\t'.join(['IV','INS',str(apos),'TBD','TBD','TBD','TBD','TBD',modified])
              out_lines.append(out_line)
            #   varNum='IV'
            else:
              # if variant is SNP or SUB
              sample_frequency = variants[spos][2]
              coverage = variants[spos][1]
              variant = w + str(apos) + r + '_' + c + '*'
              gene_region, AA_mutation = determine_codon(str(apos),c)
              out_line = '\t'.join(['IV','SNP',variant,coverage,'1',sample_frequency,gene_region,AA_mutation,modified])
              out_lines.append(out_line)
        # consensus does not equal ref or anchor
        elif c != r and r == w:
            # Type II, modified
            if c == 'n' and spos not in list(variants.keys()):
                pass
            elif r == '-':
                gene_region, AA_mutation = determine_codon(str(apos),c)
                out_line = '\t'.join(['II','INS',str(apos)+c,coverage,'1',sample_frequency,gene_region,'TBD',modified])
                out_lines.append(out_line)
            elif c == '-':
                gene_region, AA_mutation = determine_codon(str(apos),c)
                out_line = '\t'.join(['II','DEL',str(apos),coverage,'1',sample_frequency,gene_region,'TBD',modified])
                out_lines.append(out_line)
            else:
                if len(c) == 1:
                    variant = w + str(apos) + c
                    gene_region, AA_mutation = determine_codon(str(apos),c)
                    out_line = '\t'.join(['II','SNP',variant,coverage,'1',sample_frequency,gene_region,AA_mutation,modified])
                    out_lines.append(out_line)
                else:
                    print("C",c)
                    for i in range(len(c)):
                        variant = w + str(apos) + c[i]
                        gene_region, AA_mutation = determine_codon(str(apos),c[i])
                        out_line = '\t'.join(['II','SNP',variant,coverage[i],'1',sample_frequency[i],gene_region,AA_mutation,modified])
                        out_lines.append(out_line)
        # consensus does not equal anchor
        elif c == r and r != w:
            if w == '-':
                gene_region, AA_mutation = determine_codon(str(apos),c)
                out_line = '\t'.join(['I','INS',str(apos)+c,'n/a','1','1.000000',gene_region,'TBD',modified])
                out_lines.append(out_line)
            elif c == '-':
                gene_region, AA_mutation = determine_codon(str(apos),c)
                out_line = '\t'.join(['I','DEL',str(apos),'n/a','1','1.000000',gene_region,'TBD',modified])
                out_lines.append(out_line)
            # This is not a variant because both the ref and con are Ns
            elif c == 'n':
                pass
            else:
                variant = w + str(apos) + c
                gene_region, AA_mutation = determine_codon(str(apos),c)
                out_line = '\t'.join(['I','SNP',variant,'n/a','1','1.000000',gene_region,AA_mutation,modified])
                out_lines.append(out_line)
    return out_lines

# Remove trailing deletions and consolidate deletions onto one line
def process_outlines(out_lines):
    # remove opening and closing deletions when strain is shorter than anchor
    # Go backwards to remove closing deletion
    #print(f'outlines is {out_lines}')
    raw_stripped = []
    closing_deletion = True    
    # start at end of out_lines through -1 by -1
    for i in range(len(out_lines)-1,-1,-1):
        line = out_lines[i].split('\t')
        #need to redefine opening/closing deletions because if the first/final variant is an indel, it wouldn't be caught
        try:
            prev = out_lines[i+1].split('\t') #get previous entries in list so we can compare positions
            #print(f'prev is {prev}')
        except:
            prev = 'None'
            #print(f'prev is {prev}')
        #print(f'position for line should be {line} and {line[2]}')
        if prev == 'None': #skip over the first iteration
            pass
        else:
            #print(f'line in prev loop is {line}')
            #print(f'line[1] is {line[1]}, line[2] is {line[2]}, prev is {prev}')
            if closing_deletion and line[1] in ['INS','DEL'] and int(line[2]) +1 != prev[2]: #if the position of the current line is not one away from that of the previous line, we must have some region of anchor covered.
                closing_deletion = False
                raw_stripped.append(out_lines[i])
            # pass all end deletions until first variant not a deletion, then break out of loop and append variant
            elif line[1] in ['INS','DEL'] and closing_deletion:
                #print(f'process_outlines: {line} is considered a closing deletion')
                pass
            else: #handles cases where closing_deletion is already false or there is a SNP
                closing_deletion = False
                raw_stripped.append(out_lines[i])
    # GO backwards again to remove opening deletion, because backwards + backwards = forwards
    stripped = []
    opening_deletion = True
    for i in range(len(raw_stripped)-1,-1,-1):
        line = raw_stripped[i].split('\t')
        try:
            prev = raw_stripped[i+1].split('\t') #get previous entries in list so we can compare positions
        except:
            prev = 'None'
        if prev == 'None': #skip over the first iteration
            pass
        else:
            #print(f'line in prev loop is {line}')
            if opening_deletion and line[1] in ['INS','DEL'] and int(line[2]) -1 != int(prev[2]):
                opening_deletion = False
                stripped.append(raw_stripped[i])
            # pass all end deletions until first variant not a deletion, then break out of loop and append variant
            elif line[1] in ['INS','DEL'] and opening_deletion:
                pass
            else:
                opening_deletion = False
                stripped.append(raw_stripped[i])
   # for line in stripped:
       # print(line)
    # concatenate consecutive deletions and insertions not at the ends of sequence
    clean_out_lines = []
    deletion_start = ''
    deletion_end = ''
    insertion_nucs = ''
    # dict to hold co-occuring variants
    covars = {}
    # This iterates through len(stripped)-1 because it looks ahead one variant
    # to look for lines that can be concated ,e.g., several consecutive dels
    # the clean_out_lines.append(stripped[-1]) at the end of the loop add the final element
    for i in range(len(stripped)-1):
        line = stripped[i].split('\t')
        print("line: ", line)
        variant = line[1]
        if variant in ['INS','DEL'] and line[5] not in ["n/a","TBD"] and float(line[5]) < 0.5:
            clean_out_lines.append(stripped[i])
        # if the variant is a deletion
        elif variant == 'DEL':
            variant_pos = int(line[2]) # nucleotide position
            # determine where deletion begins
            if deletion_start == '':
                deletion_start = str(variant_pos)
            next_variant = stripped[i+1].split('\t')[1]
            # if the next variant is a deletion of the same type
            if next_variant == 'DEL' and line[0] == stripped[i+1].split('\t')[0]:
                next_variant_pos = int(stripped[i+1].split('\t')[2])
                # Check that the positions are in consecutive order
                if next_variant_pos == variant_pos+1:
                    # if this is a new deletion, save the end position
                    deletion_end = str(variant_pos)
                    continue
            # if deletion is 1 nt
            if deletion_end == '':
                aa_mutation = determine_del_aa_mutation(int(deletion_start),1)
                del_line = '\t'.join([line[0],'DEL',str(variant_pos),'n/a','1','1.000000',line[-3],aa_mutation,line[-1]])
                clean_out_lines.append(del_line)
                deletion_start = ''
            # if the deletion is longer than 1 nt
            else:
                start_del_end = deletion_start + '-' + str(int(deletion_end) + 1)
                del_length =  int(deletion_end) + 1 - int(deletion_start) + 1
                aa_mutation = determine_del_aa_mutation(int(deletion_start),del_length)
                del_line = '\t'.join([line[0],'DEL',start_del_end,'n/a',str(del_length),'1.000000',line[-3],aa_mutation,line[-1]])
                clean_out_lines.append(del_line)
                deletion_start = ''
                deletion_end = ''
        elif variant == 'INS': # and line[2] != "TBD"
            # line[2] == 45g
            # nuc == g
            # variant_pos = 45
            if line[2] == "TBD":
                clean_out_lines.append('\t'.join(line))
                continue
            nuc = line[2][-1]
            print("Nuc: ", nuc)
            insertion_nucs += nuc
            print("variant position: ", line[2][:-1])
            variant_pos = int(line[2][:-1]) # anchor position
            next_variant = stripped[i+1].split('\t')[1]
            # if the next variant is a insertion
            if next_variant == 'INS':
                next_variant_pos = int(stripped[i+1].split('\t')[2][:-1])
                # if the next insertion is part of the current
                # if next_variant_pos != variant_pos:
                if next_variant_pos == variant_pos: # variant position remains the same b/c is with regards to anchor positions
                    continue
            # if insertion is 1 nt long
            if len(insertion_nucs) == 1:
                aa_mutation = determine_ins_aa_mutation(insertion_nucs.upper(),variant_pos)
                insert = str(variant_pos) + '[' + nuc + ']' + str(variant_pos+1)
                ins_line = '\t'.join(['I','INS',insert,'n/a','1','1.000000',line[-3],'TBD',line[-1]])
                clean_out_lines.append(ins_line)
                insertion_nucs = ''
            # if the insert is longer than 1 nt
            else:
                aa_mutation = determine_ins_aa_mutation(insertion_nucs.upper(),variant_pos)
                insert = str(variant_pos) + '[' + insertion_nucs + ']' + str(variant_pos+1)
                ins_line = '\t'.join(['I','INS',insert,'n/a',str(len(insertion_nucs)),'1.000000',line[-3],aa_mutation,line[-1]])
                clean_out_lines.append(ins_line)
                insertion_nucs = ''
        elif variant == 'SNP':
            if stripped[i+1].split('\t')[1] != 'SNP': # check if the next variant is also a SNP, if not, pass
                pass
            else:
            # Check if two variants co-occur in same codon
                if line[0] == 'III': # if current variant is a type III
                    variant_pos = line[2][1:-8] # type III variants have _rev_nt that needs to be removed
                elif line[0] == "IV":
                    if "*" in line[2]:
                        variant_pos = line[2][1:-4]
                    else:
                        variant_pos = line[2][1:-1]
                else:
                    variant_pos = line[2][1:-1] # remove nt from front and end of variant nomenclature
                print("Variant pos: ", variant_pos)
                cur_codon = sars_genes[variant_pos]['Codon']
                cur_nuc = line[2][-1]
                cur_nuc_codon_pos = sars_genes[variant_pos]['Codon Position']
                if stripped[i+1].split('\t')[0] == 'III': # if the next SNP is also a type III variant
                    next_pos = stripped[i+1].split('\t')[2][1:-8]
                elif stripped[i+1].split('\t')[0] == 'IV':
                    if "*" in stripped[i+1].split('\t')[2]:
                        next_pos = stripped[i+1].split('\t')[2][1:-4]
                    else:
                        next_pos = stripped[i+1].split('\t')[2][1:-1]
                else:
                    next_pos = stripped[i+1].split('\t')[2][1:-1]
                print()
                print("Line Variant", line[0])
                print("Stripped: ", stripped[i+1].split('\t')[2], "Entered: ", stripped[i+1].split('\t')[0])
                print("Works 1: ", stripped[i+1].split('\t')[2][1:-8], "Works 2: ", stripped[i+1].split('\t')[2][1:-1])
                print("Variant pos: ", variant_pos, "Next pos: ", next_pos)
                print()
                next_codon = sars_genes[next_pos]['Codon']
                next_nuc = stripped[i+1].split('\t')[2][-1]
                next_nuc_codon_pos = sars_genes[next_pos]['Codon Position']
                different_position = variant_pos != next_pos
                same_gene = sars_genes[variant_pos]['Gene'] == sars_genes[next_pos]['Gene'] # sets same_gene to true or false
                same_region = sars_genes[variant_pos]['Region'] == sars_genes[next_pos]['Region']
                same_aa_num = sars_genes[variant_pos]['AA #'] == sars_genes[next_pos]['AA #']
                # if same_gene = TRUE, and same_region = TRUE, etc.
                # if 2 SNPs are in same gene, same region, same aa # and not in a UTR
                if different_position and same_gene and same_region and same_aa_num and 'UTR' not in sars_genes[variant_pos]['Gene']:
                    if int(variant_pos) > 265 or int(variant_pos) < 29675: # then it's not in the UTRs, so maybe in codon
                        # save to covars
                        aa_info = sars_genes[variant_pos]['Gene']+sars_genes[variant_pos]['Region']+sars_genes[variant_pos]['AA #']
                        if aa_info in covars.keys():
                            covars[aa_info].append(next_pos+next_nuc+next_nuc_codon_pos)
                        else:
                            covars[aa_info] = [variant_pos+cur_nuc+cur_nuc_codon_pos, next_pos+next_nuc+next_nuc_codon_pos]
            clean_out_lines.append(stripped[i])
    # This is a temporary solution that needs to be improved upon
    clean_out_lines.append(stripped[-1])
    # for line in clean_out_lines:
    #     print(line)
    # Check for SNPs that occur in the same codon and flag them
    clean_covars = {}
    # print(covars)
    # Go through covars and generate correct notation for co-occuring variants
    for k in covars.keys():
        # grab the AA and codon and initialize the mutated codon
        covar = covars[k]
        AA = sars_genes[covar[0][:-2]]['AA'] # example covar =  [212A1, 213T2] so 0 of covar = 212A1[:-2] = the aa at nt pos 212
        AAnum = sars_genes[covar[0][:-2]]['AA #']
        codon = sars_genes[covar[0][:-2]]['Codon']
        mutated_codon = [e for e in codon] # nt @ 1, nt @ 2, or nt @ 3 within anchor codon at a given nucleotide position
        # for each co-occuring variant, mutate the codon
        for var in covar:
            pos = var[:-2]
            nuc = var[-2]
            nuc_codon_pos = var[-1]
            # print(nuc,nuc_pos)
            mutated_codon[int(nuc_codon_pos)-1] = nuc
        print(mutated_codon)
        mutated_codon = ''.join(mutated_codon)
        # print(k,AA,codon,mutated_codon)
        # determine the AA of the mutated codon
        if mutated_codon.startswith("*") or "‡" in mutated_codon:
            mutated_AA = f"figure out AA for codon {mutated_codon}"
        else:
            mutated_AA = genetic_code[mutated_codon.upper()]
        #retain prev codon, prev AA, AA #, mutated AA, and mutated codon
        clean_covars[k] = [AA+AAnum+mutated_AA,codon+' -> '+mutated_codon] # Y21H, ACG -> TAG
    # print(clean_covars)
    clean_out_lines_cooccuring = []
    for i in range(len(clean_out_lines)):
        line = clean_out_lines[i].split('\t')
        variant = line[1]
        if variant == 'SNP':
            # If the SNP is in a co-occuring variant codon
            if line[0] == 'III':
                variant_pos = line[2][1:-8]
            elif line[0] == "IV":
                if "*" in line[2]:
                    variant_pos = line[2][1:-4]
                else:
                    variant_pos = line[2][1:-1]
            else:
                variant_pos = line[2][1:-1]
            aa_info = sars_genes[variant_pos]['Gene']+sars_genes[variant_pos]['Region']+sars_genes[variant_pos]['AA #']
            if aa_info in clean_covars.keys():
                new_line = list(line)
                # new_line[-2] += f" ** multiple SNPs in codon {sars_genes[variant_pos]['AA #']}: "+','.join(clean_covars[aa_info])
                new_line[-2] = f"single SNP \"{new_line[-2]}\"" + " -> combined SNP: " + f"\"{clean_covars[aa_info][0]}\""
                clean_out_lines_cooccuring.append('\t'.join(new_line))
            else:
                clean_out_lines_cooccuring.append('\t'.join(line))
        else:
            clean_out_lines_cooccuring.append('\t'.join(line))
    # clean = [clean_out_lines_cooccuring[i] for i in range(len(clean_out_lines_cooccuring)-1,-1,-1)]
    clean = list(clean_out_lines_cooccuring)
    return clean


def main():
    parser = argparse.ArgumentParser(description='Script to generate variant table \
    CSV file for report generation.')
    parser.add_argument('-f', metavar='Path to mafft alignment file.',
     type=str, nargs='?', required=True, help='The mafft alignment filepath. The order \
     of the sequences should be anchor, strain, consensus.')
    parser.add_argument('-v', metavar='Path to final VCF file',
     type=str, nargs='?', required=True)
    parser.add_argument('-o', metavar='Output tab delimited file path. Ex: output.txt',
     type=str, nargs='?', required=True, help='The desired output location file path')
    # parser.add_argument('-w', metavar='Output warning file path',
    #  type=str, nargs='?', required=False, help='The desired warning output location')
    args = parser.parse_args()
    # Check if we are using a modified strain and determine the regions that were modified
    regions = []
    if "_modified_strain" in args.o.split('/')[-2]:
        working_dir = '/'.join(args.o.split('/')[:-2])
        low_cov_n_dir = args.o.split('/')[-2].replace("_modified_strain","")
        low_cov_n_filename = low_cov_n_dir+"_low_coverage_N_regions.txt"
        print(working_dir+'/'+low_cov_n_dir+'/'+low_cov_n_filename)
        for line in open(working_dir+'/'+low_cov_n_dir+'/'+low_cov_n_filename,"r"):
            regions.append([int(e) for e in line.split('\t')])
    # Put the variants in a dictionary variants[pos] = [pos,cov,af,ref,alt]
    variants = get_variants(args.v)
    msa_path = args.f
    anchor, strain, anchor_consensus = process_alignment(msa_path)
    out_lines = generate_table(anchor,strain,anchor_consensus,variants,regions)
    out_lines = process_outlines(out_lines)
    # for line in out_lines:
    #     print(line)


    #### GO BACK THROUGH LINES AND ADD CORRECTION NOTATION FOR INDELS
    ### The alignment can generate different notation than what we have come up with
    ### in the reports
    with open(args.o,'w') as f:
        for line in out_lines:
            f.write(line+'\n')
    # with open(args.w,'w') as f:
    #     for line in warning_lines:
    #         f.write(line+'\n')


if __name__ == "__main__":
    main()
