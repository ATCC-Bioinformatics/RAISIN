import os
import sys
from dna_features_viewer import GraphicFeature, GraphicRecord
from dna_features_viewer import BiopythonTranslator
from Bio import SeqIO
from bokeh.resources import CDN
from bokeh.embed import file_html



vcf="/home/shared/service-data/863/ILM3_2023CF07_NR-59685_EPI_ISL_18507690/reads_to_reference/ILM3_2023CF07_NR-59685_reads_to_reference_EPI_ISL_18507690_final.vcf"
mafft="/home/shared/service-data/863/ILM3_2023CF07_NR-59685_EPI_ISL_18507690/reads_to_reference/ILM3_2023CF07_NR-59685_mafft_alignment.fasta"
gbk="/home/shared/references/NC_045512.2.gbk"
###! Submission 979

# vcf="/home/shared/service-data/979/ILM3_2024AL03_229E_Calu-3_Passage-2_MW202340.1/consensus_results-notds/ILM3_2024AL03_229E_Calu-3_Passage-2_reads_to_reference_MW202340.1_lofreq.vcf"

# ref="/home/shared/references/MW202340.1.fasta"
# anchor="/home/shared/references/KY684760.fasta"
# gbk="/home/shared/references/KY684760.gbk"
# mafft="/home/shared/repos/RAISIN_develop/test_mafft_alignment.fasta"

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

def get_mutated_region(region,reference,position,ref,alt):
    # set for 0 index used in biopython
    position -= 1
    # If deletion, then store deleted positions
    deletion_positions = []
    if len(ref) > len(alt):
        for i in range(1,len(ref)):
            deletion_positions.append(position+i)
    # Read genbank and pull out sequence of region adding in the variant
    for r in SeqIO.parse(open(gbk,"r"), "genbank"):
        for f in r.features:
            if f.type == 'CDS':
                try:
                    label = f.qualifiers['product']
                except:
                    label = f.qualifiers['note']
                print(label)
                # print("Supposed product: ", f.qualifiers['product'])
#            if f.type == 'CDS' and f.qualifiers['product'][0] == region:
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
                        # if l < 40:
                        #     print(l, position, "Ref: ", ref, "Alt: ", alt)
                        #     print(nt_seq)
                    if f.location.strand == -1:
                        nt_seq = ''.join([revcom[e] for e in nt_seq[::-1]])
                    break
    return nt_seq


def determine_variant_type(ref, alt, cons):
    if ref != alt != cons:
        var_type = "IV"
    elif ref == alt != cons:
        var_type = "II"
    elif ref == cons != alt:
        var_type = "III"
    elif ref != alt == cons:
        var_type = "I"
    return var_type



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

anchor, strain, sample_consensus = process_alignment(mafft)

### Anchor Mode
variants = {}
for line in open(vcf):
    if "#" not in line:
        line = line.strip().split("\t")
        pos = line[1]
        ref = line[3]
        alt = line[4]
        af = [e.replace("AF=","") for e in line[7].split(';') if "AF" in e][0]
        if pos in variants.keys():
            num_alts +=1
            pos = pos + "_"+str(num_alts)
            variants[pos] = [ref,alt,af]
        else:
            num_alts = 1
            variants[pos] = [ref,alt,af]

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

## if the vcf position
# {0: [-1], 1: [-1], 2: [-1], 3: [-1], 4: 0, 5: 1, 6: 2, 7: 3, 8: 4, 9: [4], 10: [4], 11: [4], 12: 5, 13: 6, '13-14': [7, 8, 9], 14: 10}
anchor_strain_dict = anchor_genomes(anchor, strain)
anchor_consen_dict = anchor_genomes(anchor, sample_consensus)

vcf_positions = list(variants.keys())
final_variants = {k:variants[k] for k in variants.keys()}
for i in range(len(vcf_positions)):
    vcf_pos = vcf_positions[i]
    for anchor_pos, strain_pos in anchor_strain_dict.items():
        if type(strain_pos) is not list: ## address non-indels first (positions of indels are in a list format)
            if int(vcf_pos) == int(strain_pos):
                print(vcf_pos, strain_pos, anchor_pos)
                final_variants[vcf_pos].append(anchor_pos)
        if type(anchor_pos) == str and anchor_pos.count("-") == 1:
            for nuc in strain_pos:
                if int(vcf_pos) == int(nuc):
                    anchor_start = int(anchor_pos.split("-")[0])+1
                    final_variants[vcf_pos].append(anchor_pos)

reference = ""
graphic_features = []
features = {}
for r in SeqIO.parse(open(gbk), "genbank"):
    if reference == "":
        reference = ''.join(list(r.seq))
    for f in r.features:
        # print("Type: ", f.type)
        if f.type == "source":
            org_id = gbk.replace(".gb","")
            if 'organism' in f.qualifiers:
                org_id += f" {f.qualifiers['organism'][0]}"
            if 'strain' in f.qualifiers:
                org_id += f" {f.qualifiers['strain'][0]}"
        elif f.type == 'CDS':
            # print(f.location)
            # print(f.location.strand)
            # for splice sites, the location shows up as join(15..44,517..852) instead of 15..674, need to account for cases where variant is between splice site, i.e. 122
            for i in range(0,len(f.location.parts)):  
                try:
                    label = f.qualifiers['product'][0]
                except:
                    label = f.qualifiers['note'][0]
                label = label + str(i)  # nuclear export protein 0, nuclear export protein 1
                start = str(f.location.parts[i]).split('[')[1].split(']')[0].split(':')[0]
                end = str(f.location.parts[i]).split('[')[1].split(']')[0].split(':')[1]
                # print("Label:", label, start, end)
                features[label] = {}
                features[label]['start'] = start
                features[label]['end'] = end
                features[label]['reference2regionposition'] = {}
                features[label]['strand'] = f.location.strand
                aa_seq = f.qualifiers['translation'][0]
                nt_seq = ""
                if f.location.strand == 1:
                    c=0
                    for l in f.location:
                        nt_seq += r.seq[l]
                        features[label]['reference2regionposition'][l] = c
                        c+=1
                elif f.location.strand == -1:
                    c=0
                    for l in f.location:
                        nt_seq += r.seq[l]
                        features[label]['reference2regionposition'][l] = c
                        c+=1
                    nt_seq = ''.join([revcom[e] for e in nt_seq])
                # don't include stop codon
                features[label]['aa_seq'] = aa_seq
                features[label]['nt_seq'] = nt_seq
            ### For ensuring correct mucleotide sequence and translation
            # aas = ''.join([genetic_code(seq[i:i+3].upper()) for i in range(0,len(seq),3)])
            # if aas[:-1] == f.qualifiers['translation'][0]:
                # print("correct translation")

# {0: [-1], 1: [-1], 2: [-1], 3: [-1], 4: 0, 5: 1, 6: 2, 7: 3, 8: 4, 9: [4], 10: [4], 11: [4], 12: 5, 13: 6, '13-14': [7, 8, 9], 14: 10}

anchor_seq=anchor.replace("-", "")
strain_seq=strain.replace("-", "")
consen_seq=sample_consensus.replace("-", "")

##! TO DO: incorporate UTR variants found when doing an MSA
# Grab the start position of where the first protein starts (to know where 5'UTR ends)
first_start = min([int(features[k]['start']) for k in features.keys()])
# Grab the end position of where the last protein ends (to know where 3'UTR begins)
last_end = max([int(features[k]['end']) for k in features.keys()])

deletion_pos = []
msa_variants = {}
for pos in anchor_strain_dict.keys():
    for k in features.keys(): #where features is a parsed gbk)
        prot_start = int(features[k]['start'])
        prot_end = int(features[k]['end'])
        alt_pos = anchor_strain_dict[pos]
        cons_pos = anchor_consen_dict[pos]
        region = features[k]['nt_seq']
        region_translation = features[k]['aa_seq']
        protein = k[:-1]
        ##* INSERTION!! # insertions are represented like this in the dict: {4-5: [1, 2, 3, 4]}, we want to change it to strain_seq[1:5]
        if type(pos) == str and "-" in pos: # if position is a string, i.e. "-", then we're in an insertion, i.e. "3-4" or "-1-0"
            print("Insertion")
            if pos.count("-") == 2: #.e. "-1-0"
                print("WARNING: The reference genome hasn't started yet. Looks like the strain reference might be longer than the anchor reference.")
                break
            else:
                ref_start = int(pos.split("-")[0])
                ref = anchor_seq[ref_start]
            if ref_start > prot_start and ref_start < prot_end: # check region based on ref position
                one_before_ins = alt_pos[0] - 1 # nucleotide right before start of insertion
                ins_start = alt_pos[0] + 1 # actual position of insertion start (1-based)
                ins_end = alt_pos[-1] + 1 # actual position of insertion end (1-based)
                alt = strain_seq[one_before_ins:ins_end]
                cons = consen_seq[one_before_ins:ins_end]
                insertion_nucleotides = strain_seq[alt_pos[0]:ins_end]
                variant_type = determine_variant_type(ref, alt, cons)
                msa_variants[ref_start+1] = ["INS", variant_type, ref, alt, cons, "1.000000"]
                mutated_region = get_mutated_region(protein,reference,ref_start,ref,alt)
                mutated_translation = ""
                for i in range(0,len(mutated_region),3):
                    aa = genetic_code(mutated_region[i:i+3].upper())
                    if i+3>len(mutated_region)-1 or aa == "*":
                        break
                    else:
                        mutated_translation+=aa
                msa_variants[ref_start+1].append(protein)
                msa_variants[ref_start+1].append(indel_notation(region_translation,mutated_translation,len(insertion_nucleotides)%3!=0,"ins"))
        elif alt_pos == [-1]:
            print("Strain reference hasn't started, skipping")
            break
        elif pos > prot_start and pos < prot_end:
            ## * DELETION!! # {8:4, 9: [4], 10: [4], 11: [4], 12:5}
            if type(alt_pos) == list:
                print("Deletion")
                # Grab all the positions of the deletion
                deletion_pos.append(pos)
                if pos+1 not in anchor_strain_dict.keys():
                    # print("Position at end of reference!")
                    break
                elif alt_pos == anchor_strain_dict[pos+1]: #if we're still in the deletion, keep moving through
                    continue
                else: # we've reached the end of the deletion
                    one_before_del = deletion_pos[0] - 1 # nucleotide right before start of deletion
                    del_start = deletion_pos[0] + 1 # actual position of deletion start (1-based)
                    del_end = deletion_pos[-1] + 1 # actual position of deletion end (1-based)
                    deleted_nucs = anchor_seq[deletion_pos[0]:del_end]
                    ref = anchor_seq[one_before_del:del_end]
                    alt = strain_seq[alt_pos[0]]
                    cons = strain_seq[cons_pos[0]]
                    variant_type = determine_variant_type(ref, alt, cons)
                    mutated_region = get_mutated_region(protein,reference,del_start,ref,alt)
                    mutated_translation = ""
                    msa_variants[del_start] = ["DEL", variant_type, ref, alt, cons, "1.000000"]
                    for i in range(0,len(mutated_region),3):
                        aa = genetic_code(mutated_region[i:i+3].upper())
                        if i+3>len(mutated_region)-1 or aa == "*":
                            break
                        else:
                            mutated_translation+=aa
                    msa_variants[del_start].append(protein)
                    msa_variants[del_start].append(indel_notation(region_translation,mutated_translation,len(deleted_nucs)%3!=0,"del"))
                    deletion_pos = [] ## reset deletion pos
                    break
            else:
                ref = anchor_seq[pos]
                alt = strain_seq[alt_pos]
                cons = consen_seq[cons_pos]
                full_position = pos + 1  # actual position of SNP (1-based)
                if ref == alt == cons: # if there is no snp
                    break
                if alt.lower() == "n":
                    break
                variant_type = determine_variant_type(ref, alt, cons)
                if variant_type == "II": # if consensus is the nucleotide of difference
                    mutated_region = get_mutated_region(protein,reference,pos+1,ref,cons)
                else:
                    mutated_region = get_mutated_region(protein,reference,pos+1,ref,alt)
                if mutated_region != region:
                    codon_number = 0
                    for i in range(0,len(mutated_region),3):
                        codon_number+=1
                        if mutated_region[i:i+3] != region[i:i+3]:
                            original_codon = region[i:i+3]
                            mutated_codon = mutated_region[i:i+3]
                            original_aa = genetic_code(original_codon.upper())
                            mutated_aa = genetic_code(mutated_codon.upper())
                            break
                            # print(c,mutated_region[i:i+3], region[i:i+3])
                    msa_variants[full_position] = ["SNP", variant_type, ref, alt, cons, "1.000000"]
                    # add notation to VCF
                    if original_aa == mutated_aa:
                        msa_variants[full_position].append(protein)
                        msa_variants[full_position].append("--")
                    else:
                        msa_variants[full_position].append(protein)
                        msa_variants[full_position].append(original_aa+str(codon_number)+mutated_aa)
        # else:
        #     print("Something weird is going on")
        #     print(pos, protein, "start: ", int(features[k]['start']), "end: ", int(features[k]['end']))
        # break

#! Merge with vcf to get frequency information
tv = {k:msa_variants[k] for k in msa_variants.keys()}
for anchor_pos in tv.keys():
    for var in final_variants.keys():
        anchor_var = final_variants[var][-1]
        if anchor_pos == anchor_var:
            print("Anchor pos: ", anchor_var, tv[anchor_pos])
            frequency = final_variants[var][2]
            cons_allele = final_variants[var][1]
            tv[anchor_pos][5] = frequency
            tv[anchor_pos][4] = cons_allele


with open(f"tv_tester_variants.txt","w") as f:
    f.write("\t".join(["Anchor Position","Variant", "Variant Type", "Anchor Allele","Strain Allele", "Sample Allele",
        "Allele Frequency", "Protein","AA mutation"]))
    f.write("\n")
    for k in tv.keys():
        f.write(str(k)+"\t"+"\t".join(tv[k])+"\n")



with open(f"tester_variants.txt","w") as f:
    f.write("\t".join(["Anchor Position","Variant", "Variant Type", "Anchor Allele","Strain Allele", "Sample Allele",
        "Allele Frequency", "Protein","AA mutation"]))
    f.write("\n")
    for k in msa_variants.keys():
        f.write(str(k)+"\t"+"\t".join(msa_variants[k])+"\n")

#     ###! Type I indel - all nucs in strain and consensus match, diff to anchor
#     ###! Type II indel - strain and consensus differ, but strain and anchor agree
#     ###! Type III indel - all nucs in anchor and consensus match, diff to strain
#     ###! Type IV indel - strain, consensus, and anchor differ


###########
"""
test_d = {k: anchor_strain_dict[k] for k in range(11284, 11300)}
deletion_pos = []
for k in test_d.keys():
    alt_pos = test_d[k]
    print("On position", k, "Actual position: ", k+1, alt_pos)
    if type(test_d[k]) == list:
        print("Deletion")
        print("Anchor: ", anchor[k], "Strain: ", strain[k])
        deletion_pos.append(k)
        if k+1 not in test_d.keys():
            print("At end of reference")
        elif alt_pos == test_d[int(k)+1]:
            continue
        else:
            print("Saving deletion at pos ", k)
            one_before_del = deletion_pos[0] - 1
            del_start = deletion_pos[0] + 1
            del_end = deletion_pos[-1] + 1
            print("Start: ", del_start, "End: ", del_end)
            deleted_nucs = anchor_seq[deletion_pos[0]:del_end]
            ref = anchor_seq[one_before_del:del_end]
            alt = strain_seq[alt_pos[0]]
            print("Ref: ", ref, "Alt: ", alt)
            deletion_pos = []
    else:
            print("Anchor: ", anchor[k], "Strain: ", strain[k])
    print()

"""