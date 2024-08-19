### Functions used in universal_raisin.py

import os
import sys
from Bio import SeqIO

revcom = {"A":"T","T":"A","G":"C","C":"G"}
#######################* Function Time!! 
## Genetic Code ----------------------------------------
def genetic_code(codon):
    table = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
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
    try:
        return table[codon]
    except:
        return 'X'

def replace_degenerate_nucleotides(consensus_nuc, strain_nuc):
    degen_table = {
        'W': ['A', 'T'], 'S': ['C', 'G'], 'R': ['A', 'G'],
        'Y': ['C', 'T'], 'K': ['G', 'T'], 'M': ['A', 'C'],
    }
    try:
        consen = consensus_nuc.upper()
        strain = strain_nuc.upper()
        choices = degen_table[consen]
        # Since consensus nucleotide is a degenerate nuc, the true identity
        # of the nucleotide should be the nuc that does not match the strain
        for nuc in choices:
            if strain != nuc:
                return nuc
    except:
        return consensus_nuc

def get_mutated_region(region,reference,position,ref,alt, gbk_file):
    # set for 0 index used in biopython
    position -= 1
    # If deletion, then store deleted positions
    deletion_positions = []
    if len(ref) > len(alt):
        for i in range(1,len(ref)):
            deletion_positions.append(position+i)
    # Read genbank and pull out sequence of region adding in the variant
    for r in SeqIO.parse(open(gbk_file,"r"), "genbank"):
        for f in r.features:
            if f.type == 'CDS':
                try:
                    label = f.qualifiers['product']
                except:
                    label = f.qualifiers['note']
                if label[0] == region:
                    nt_seq = ""
                    # print(f)
                    for l in sorted(list(f.location)):
                        # If the position is a deleted position, don't add nt
                        if l in deletion_positions:
                            pass
                        elif l == position:
                            # If the variant is a SNP, add the variant
                            if len(ref) == 1 and len(alt) == 1:
                                nt_seq += alt 
                            # If the variant is an INS, add the entire alt
                            elif len(alt) > len(ref):
                                nt_seq += alt 
                        # Otherwise, add the reference nt
                        else:
                            nt_seq += reference[l]
                    if f.location.strand == -1:
                        nt_seq = ''.join([revcom[e] for e in nt_seq[::-1]])
                    break
    return nt_seq

def indel_notation(original,mutated,frameshift,var):
    if frameshift == False:
        original_diff_aas = original 
        mutated_diff_aas = mutated 
        # Increment through each sequence, if the nts are the same, then remove the
        # first position from the *diff* sequence, this leaves the parts of the 
        # original and mutated sequences from the first different nt through
        # the end of the sequences, 
        # e.g., ATCTAC and ATCGAC -> TAC GAC (first three removed because they are the same)
        for i in range(len(original)):
            if i > len(mutated)-1:
                break
            elif original[i] == mutated[i]:
                original_diff_aas = original_diff_aas[1:]
                mutated_diff_aas = mutated_diff_aas[1:]
            else:
                break
        # Increment backwards through the sequences removing nts that are the same
        # similar to above. 
        # e.g., TAC GAC -> T G (last two nts removed)
        # Only do this if both *diff* sequences contain nts, otherwise the mutation
        # generated an early stop codon or removed a stop codon
        if len(mutated_diff_aas) > 0 and len(original_diff_aas) > 0:
            for j in range(1,len(original)):
                if original[-j] != mutated[-j]:
                    if j > 1:
                        original_diff_aas = original_diff_aas[:-j+1]
                        mutated_diff_aas = mutated_diff_aas[:-j+1] 
                    break
        if var == "del":
            if len(mutated_diff_aas) > 0 and len(original_diff_aas) > 0:
                notation = "Δ" + original_diff_aas[:-1] + f" (aa{i+1}-{i+len(original_diff_aas[:-1])}), " 
                notation += original_diff_aas[-1] + str(j) + mutated_diff_aas
                return notation
            # mutation caused insetion of early stop codon
            elif len(original_diff_aas) > 0:
                notation = "Δ" + original_diff_aas + f" (aa{i+1}-{i+len(original_diff_aas)})" 
                return notation                 
            # mutation removed stop codon and extended protein
            else:
                if len(mutated_diff_aas) > 1:
                    notation = "*"+str(len(original)+1)+ "["+mutated_diff_aas +"*]"
                else:
                    notation = "*"+str(len(original)+1)+mutated_diff_aas[0]
                return notation 
        elif var == "ins":
            notation = f"{i}[" + mutated_diff_aas + f"]{i+1}" 
            return notation             
    else:
        for i in range(len(original)):
            if original[i] != mutated[i] or i+1 > len(mutated)-1 or i+1 > len(original)-1:
                break
        # If the mutation is very long, don't include the amino acids in variants.txt
        if len(original[i:]) > 10 or len(mutated[i:]) > 10:
            notation = f"aa{i+1}-{i+len(original[i:])} " + " -> "
            notation += f"aa{i+1}-{i+len(mutated[i:])} " + mutated[i:]
        else:
            notation = f"aa{i+1}-{i+len(original[i:])} " + original[i:] + " -> "
            notation += f"aa{i+1}-{i+len(mutated[i:])} " + mutated[i:]
        return notation

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

def determine_variant_type(ref, alt, cons):
    if ref != alt != cons and ref != cons:
        var_type = "IV"
    elif ref == alt != cons:
        var_type = "II"
    elif ref == cons != alt:
        var_type = "III"
    elif ref != alt == cons:
        var_type = "I"
    return var_type

def anchor_genomes(anchor_sequence, sequence):
    apos = 0        # anchor position
    spos = 0        # strain position
    position_dict = {}
    spos_list = []
    for i in range(len(anchor_sequence)):
        # NTs
        a_nuc = anchor_sequence[i]
        s_nuc = sequence[i]
        # this block specifically deals with insertions in strain/deletion in reference
        if a_nuc == "-":
            spos_list.append(spos) # if a_nuc is a dash, add the position of the inserted nucleotide in strain
            # we need to someway to tell that the insertion is done (that's when we'll add it to the dictionary)
            # check the next nucleotide, if that is not a dash, then the insertion is over. Exception if the indel is at the end of genome
            if i+1 >= len(anchor_sequence): # if the next position does not exist
                position_dict[apos] = spos_list
                spos_list = []
            else: # if the next nucleotide does exist, we're not at the end of the sequence
                a_nuc_next = anchor_sequence[i+1]
                if a_nuc_next != '-':
                    # since this position does not technically exist in the reference, the key will be the positions
                    # that the insert is in between, i.e 4-5
                    insert_pos = f"{apos-1}-{apos}"
                    position_dict[insert_pos] = spos_list
                    spos_list = [] # reset spos_lis  
        # Otherwise, proceed as normal for the anchor reference 
        if a_nuc != '-':
            # One expection, when we encounter a deletion in the strain
            # For this, we will record the position of the nucleotide before
            if s_nuc == '-':
                position_dict[apos] = [spos-1]
            else:
                position_dict[apos] = spos
            apos+=1
        if s_nuc != '-':
            spos+=1
    return position_dict

