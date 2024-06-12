## INPUTS: VCF, GBK, MODE, output_name, and output_dir, and mafft (ANCHOR MODE ONLY)
## Output: amino acid mutations

import os
import sys
from Bio import SeqIO
from universal_functions import *

# Import in sys arguments passed from run_raisin
vcf_file = sys.argv[1]
gbk_file = sys.argv[2]
mode = sys.argv[3]
output_name = sys.argv[4]
output_dir = sys.argv[5]

if len(sys.argv) == 7: # only applicable for anchor samples
    mafft = sys.argv[6]


#######* Parse GBK ##############################
revcom = {"A":"T","T":"A","G":"C","C":"G"}
reference = ""
graphic_features = []
features = {}
for r in SeqIO.parse(open(gbk_file,"r"), "genbank"):
    if reference == "":
        reference = ''.join(list(r.seq))
    for f in r.features:
        if f.type == "source":
            org_id = gbk_file.replace(".gb","")
            if 'organism' in f.qualifiers:
                org_id += f" {f.qualifiers['organism'][0]}"
            if 'strain' in f.qualifiers:
                org_id += f" {f.qualifiers['strain'][0]}"
        elif f.type == 'CDS':
            # for splice sites, the location shows up as join(15..44,517..852) instead of 15..674, need to account for cases where variant is between splice site, i.e. 122
            for i in range(0,len(f.location.parts)):  
                try:
                    label = f.qualifiers['product'][0]
                except:
                    label = f.qualifiers['note'][0]
                label = label + str(i)  # nuclear export protein 0, nuclear export protein 1
                start = str(f.location.parts[i]).split('[')[1].split(']')[0].split(':')[0]
                end = str(f.location.parts[i]).split('[')[1].split(']')[0].split(':')[1]
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

##############################################################################*

#######* Parse VCF ##############################
variants = {}
for line in open(vcf_file):
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

## If there are no variants, exit since there is nothing else to do.
# Only applies to STANDARD mode (as in ANCHOR, we can pull MSA variants)
if len(list(variants.keys())) == 0 and mode == "STANDARD":
    print("No Variants in VCF")
    quit()

# convert vcf strain position to anchor positions
if mode == "ANCHOR":
    anchor, strain, sample_consensus = process_alignment(mafft)
    anchor_strain_dict = anchor_genomes(anchor, strain)
    anchor_consen_dict = anchor_genomes(anchor, sample_consensus)
    anchor_variants = {}
    strain_positions = list(variants.keys())
    for i in range(0, len(strain_positions)):
        vcf_pos = strain_positions[i] # each position in vcf
        ref_nuc = variants[vcf_pos][0]
        alt_nuc = variants[vcf_pos][1]
        print("Position: ", vcf_pos, ref_nuc, alt_nuc)
        print(len(ref_nuc), len(alt_nuc))
        for anchor_pos, strain_pos in anchor_strain_dict.items():
                if type(strain_pos) is not list:
                    if int(vcf_pos) == int(strain_pos):
                        if len(ref_nuc) == 1: ## address SNPs and insertion first 
                            anchor_variants[str(anchor_pos)] = variants[vcf_pos]
                        elif len(ref_nuc) > 1: ##Deletions
                            anchor_variants[str(anchor_pos+1)] = variants[vcf_pos]
    variants = anchor_variants

##############################################################################*


#***************************** Variant processing!
#* STANDARD VCF Processing (Do this for ANCHOR as well)
# Go through variants
for position in variants.keys():
    full_position = position
    # "_" is present if there are multiple variants at a position
    if "_" in position:
        position = position[:position.index("_")]
    if int(position) == -10:
        print("This should never happen. It's just a chunk of code to keep things working.")
    else:
        for k in features.keys():
            start = int(features[k]['start'])
            end = int(features[k]['end'])
            protein = k[:-1] # since we added the part number to the label earlier, remove it now so it can match correctly 
            # print(start, end, k, protein)
            # If it is, then we know the protein/gene and we need to determine the aa mutation(s)
            if int(position) >= start and int(position) <= end:
                # To get the aa mutation(s) pull out the reference CDS and generate the consensus CDS
                region = features[k]['nt_seq']
                # If the lengths of the ref and alt alleles are the same and == 1 then SNP
                if len(variants[position][0])==len(variants[position][1]) and len(variants[position][1]) == 1:
                    # print(protein,reference,int(position),variants[position][0],variants[position][1])
                    mutated_region = get_mutated_region(protein,reference,int(position),variants[position][0],variants[position][1],gbk_file)
                    if mutated_region == region:
                        print("mutated and region are same")
                        variants[full_position].append('UTR')
                        variants[full_position].append("n/a")
                    else:
                        codon_number = 0
                        for i in range(0,len(mutated_region),3):
                            codon_number+=1
                            # print(i, i+3, mutated_region[i:i+3], region[i:i+3])
                            if mutated_region[i:i+3] != region[i:i+3]:
                                original_codon = region[i:i+3]
                                mutated_codon = mutated_region[i:i+3]
                                original_aa = genetic_code(original_codon.upper())
                                mutated_aa = genetic_code(mutated_codon.upper())
                                break
                                # print(c,mutated_region[i:i+3], region[i:i+3])
                        # add notation to VCF
                        if original_aa == mutated_aa:
                            variants[full_position].append(protein)
                            variants[full_position].append("--")
                        else:
                            variants[full_position].append(protein)
                            variants[full_position].append(original_aa+str(codon_number)+mutated_aa)
                # If ref allele length > alt allele length then DEL
                elif len(variants[position][0])>len(variants[position][1]):
                    # do not include stop codon
                    region_translation = features[k]['aa_seq']
                    # mutate the reference to include the deletion
                    deletion_length = len(variants[position][0])-len(variants[position][1])
                    mutated_region = get_mutated_region(protein,reference,int(position),variants[position][0],variants[position][1],gbk_file)
                    mutated_translation = ""
                    for i in range(0,len(mutated_region),3):
                        aa = genetic_code(mutated_region[i:i+3].upper())
                        if i+3>len(mutated_region)-1 or aa == "*":
                            break
                        else:
                            mutated_translation+=aa
                    variants[full_position].append(protein)
                    variants[full_position].append(indel_notation(region_translation,mutated_translation,deletion_length%3!=0,"del"))
                # If ref allele length > alt allele length then INS
                elif len(variants[position][0])<len(variants[position][1]):
                    # do not include stop codon
                    region_translation = features[k]['aa_seq']
                    # mutate the reference to include the insertion
                    insertion_nucleotides = variants[position][1]
                    mutated_region = get_mutated_region(protein,reference,int(position),variants[position][0],variants[position][1],gbk_file)
                    mutated_translation = ""
                    for i in range(0,len(mutated_region),3):
                        aa = genetic_code(mutated_region[i:i+3].upper())
                        if i+3>len(mutated_region)-1 or aa == "*":
                            break
                        else:
                            mutated_translation+=aa
                    variants[full_position].append(protein)
                    variants[full_position].append(indel_notation(region_translation,mutated_translation,len(insertion_nucleotides)%3!=0,"ins"))
                break
    if len(variants[full_position])<5:
        variants[full_position].append("UTR")
        variants[full_position].append("n/a")

######## Write to variants.txt ##########
with open(f"{output_dir}/{output_name}_variants.txt","w") as f:
    f.write("\t".join(["Position","Reference allele","Alternate allele","Allele frequency","Protein","AA mutation"]))
    f.write("\n")
    for k in variants.keys():
        f.write(k+"\t"+"\t".join(variants[k])+"\n")


######## Write to low_coverage.txt ##########
with open(f"{output_dir}/{output_name}_lowcoverage.txt","w") as f:
    f.write("\t".join(["Chromosome","Start_0index","End_1index","Average_coverage","Protein"]))
    f.write("\n")
    bed_path = sys.argv[1].split("/")
    bed_path = "/".join(bed_path[:-1])+"/low-coverage-regions.bed"
    for line in open(bed_path,"r"):
        start = line.split("\t")[1]
        end = line.split("\t")[2]
        region = []
        for k in features.keys():
            fstart = int(features[k]['start'])
            fend = int(features[k]['end'])
            protein = k[:-1]
            # If it is, then we know the protein/gene and we need to determine the aa mutation(s)
            if (int(start) >= fstart and int(start) <= fend) or (int(end) >= fstart and int(end) <= fend):
                region.append(protein)
                
        if len(region) > 0:
            f.write(line.strip()+f"\t{';'.join(region)}\n")
        else:
            f.write(line.strip()+"\tUTR\n")

## Find MSA variants
if mode == "ANCHOR":

    # Removing the dashes allows us to parse through the unaligned sequences    
    anchor_seq=anchor.replace("-", "")
    strain_seq=strain.replace("-", "")
    consen_seq=sample_consensus.replace("-", "")

    # Since the gbk don't have UTR regions, we need to put them in manually
    # To do so, figure out where the first region begins and the last region ends
    first_start = min([int(features[k]['start']) for k in features.keys()])
    last_end = max([int(features[k]['end']) for k in features.keys()])
    # all regions have a trailing 0 when the gbk is parsed, so mimic that
    features["5'UTR0"] = {'start': '0', 'end': str(first_start-1), 
        'reference2regionposition': '', 'strand': '', 'aa_seq': 'n/a', 'nt_seq': "n/a"}
    features["3'UTR0"] = {'start': str(last_end+1), 'end': str(len(reference)-1), 
        'reference2regionposition': '', 'strand': '', 'aa_seq': 'n/a', 'nt_seq': "n/a"}

    deletion_pos = []
    strain_deletion = []
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
                    if "'UTR" in protein:
                        msa_variants[ref_start+1] = ["INS", variant_type, ref.upper(), alt.upper(), cons.upper(), "1.000000", protein, "n/a"]
                    else:
                        msa_variants[ref_start+1] = ["INS", variant_type, ref.upper(), alt.upper(), cons.upper(), "1.000000"]
                        mutated_region = get_mutated_region(protein,reference,ref_start,ref,alt,gbk_file)
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
                    print("Deletion time")
                    print("Anchor pos: ", pos, "Alt pos: ", alt_pos, "Consen pos: ", cons_pos)
                    # Grab all the positions of the deletion
                    deletion_pos.append(pos)
                    if pos+1 not in anchor_strain_dict.keys():
                        # print("Position at end of reference!")
                        break
                    elif alt_pos == anchor_strain_dict[pos+1]: #if we're still in the deletion, keep moving through
                        print("Continuing", alt_pos)
                        continue
                    else: # we've reached the end of the deletion
                        one_before_del = deletion_pos[0] - 1 # nucleotide right before start of deletion
                        del_start = deletion_pos[0] + 1 # actual position of deletion start (1-based)
                        del_end = deletion_pos[-1] + 1 # actual position of deletion end (1-based)
                        deleted_nucs = anchor_seq[deletion_pos[0]:del_end]
                        ref = anchor_seq[one_before_del:del_end]
                        alt = strain_seq[alt_pos[0]]
                        cons = consen_seq[cons_pos[0]]
                        variant_type = determine_variant_type(ref, alt, cons)
                        if "'UTR" in protein:
                            msa_variants[del_start] = ["DEL", variant_type, ref.upper(), alt.upper(), cons.upper(), "1.000000", protein, "n/a"]
                        else:
                            mutated_region = get_mutated_region(protein,reference,del_start,ref,alt,gbk_file)
                            mutated_translation = ""
                            msa_variants[del_start] = ["DEL", variant_type, ref.upper(), alt.upper(), cons.upper(), "1.000000"]
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
                elif type(cons_pos) == list and type(alt_pos) != list:
                    print("Type II Deletion")
                    # Grab all the positions of the deletion
                    deletion_pos.append(pos)
                    strain_deletion.append(alt_pos)
                    if pos+1 not in anchor_consen_dict.keys():
                        # print("Position at end of reference!")
                        break
                    elif cons_pos == anchor_consen_dict[pos+1]: #if we're still in the deletion, keep moving through
                        continue
                    else: # we've reached the end of the deletion
                        one_before_del = deletion_pos[0] - 1 # nucleotide right before start of deletion
                        del_start = deletion_pos[0] + 1 # actual position of deletion start (1-based)
                        del_end = deletion_pos[-1] + 1 # actual position of deletion end (1-based)
                        deleted_nucs = anchor_seq[deletion_pos[0]:del_end]
                        ref = anchor_seq[one_before_del:del_end]
                        ## Do the same to get positions for the strain
                        strain_one_before_del = strain_deletion[0] - 1 # nucleotide right before start of deletion
                        strain_del_start = strain_deletion[0] + 1 # actual position of deletion start (1-based)
                        strain_del_end = strain_deletion[-1] + 1 # actual position of deletion end (1-based)
                        alt = strain_seq[strain_one_before_del:strain_del_end]
                        cons = consen_seq[cons_pos[0]]
                        variant_type = determine_variant_type(ref, alt, cons)
                        if "'UTR" in protein:
                            msa_variants[del_start] = ["DEL", variant_type, ref.upper(), alt.upper(), cons.upper(), "1.000000", protein, "n/a"]
                        else:
                            mutated_region = get_mutated_region(protein,reference,del_start,ref,cons,gbk_file)
                            mutated_translation = ""
                            msa_variants[del_start] = ["DEL", variant_type, ref.upper(), alt.upper(), cons.upper(), "1.000000"]
                            for i in range(0,len(mutated_region),3):
                                aa = genetic_code(mutated_region[i:i+3].upper())
                                if i+3>len(mutated_region)-1 or aa == "*":
                                    break
                                else:
                                    mutated_translation+=aa
                            msa_variants[del_start].append(protein)
                            msa_variants[del_start].append(indel_notation(region_translation,mutated_translation,len(deleted_nucs)%3!=0,"del"))
                            deletion_pos = [] ## reset deletion pos
                            strain_deletion = []
                            break
                else:
                    print("Anchor pos: ", pos, "Alt pos: ", alt_pos, "Consen pos: ", cons_pos)
                    ref = anchor_seq[pos]
                    alt = strain_seq[alt_pos]
                    cons = consen_seq[cons_pos]
                    cons = replace_degenerate_nucleotides(cons, alt)
                    full_position = pos + 1  # actual position of SNP (1-based)
                    if ref == alt == cons: # if there is no snp
                        break
                    if alt.lower() == "n":
                        break
                    if cons.lower() == "n":
                        break
                    variant_type = determine_variant_type(ref, alt, cons)
                    if "'UTR" in protein:
                        msa_variants[full_position] = ["SNP", variant_type, ref.upper(), alt.upper(), cons.upper(), "1.000000", protein, "n/a"]
                    else:
                        if variant_type == "II": # if consensus is the nucleotide of difference
                            mutated_region = get_mutated_region(protein,reference,pos+1,ref,cons,gbk_file)
                        else:
                            mutated_region = get_mutated_region(protein,reference,pos+1,ref,alt,gbk_file)
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
                            msa_variants[full_position] = ["SNP", variant_type, ref.upper(), alt.upper(), cons.upper(), "1.000000"]
                            # add notation to VCF
                            if original_aa == mutated_aa:
                                msa_variants[full_position].append(protein)
                                msa_variants[full_position].append("--")
                            else:
                                msa_variants[full_position].append(protein)
                                msa_variants[full_position].append(original_aa+str(codon_number)+mutated_aa)

    #! Merge with vcf to get frequency information for type II variants
    for anchor_pos in msa_variants.keys():
        for anchor_var in variants.keys():
            if int(anchor_pos) == int(anchor_var):
                # print("Anchor pos: ", anchor_var, msa_variants[anchor_pos])
                frequency = variants[anchor_var][2]
                cons_allele = variants[anchor_var][1]
                msa_variants[anchor_pos][5] = frequency
                msa_variants[anchor_pos][4] = cons_allele
    
    # We need to add any missing vcf variants to the msa
    msa_anchors = list(msa_variants.keys())
    var_anchors = list(variants.keys())
    # Anchor positions found in vcf that wasn't found by the msa
    not_in_msa = [x for x in var_anchors if int(x) not in msa_anchors]
    for anchor_position in not_in_msa:
        # Get the anchor position and its corresponding nucleotide
        anchor_pos = int(anchor_position)
        anchor_nuc = anchor_seq[anchor_pos-1].upper()
        # Grab everything we need from the parsed vcf
        ref = variants[anchor_position][0].upper()
        cons = variants[anchor_position][1].upper()
        freq = variants[anchor_position][2]
        region = variants[anchor_position][3]
        aa_mut = variants[anchor_position][4]
        # Determine the variant type, it should either be Type II or Type III
        variant_type = determine_variant_type(anchor_nuc, ref, cons)
        # Determine the type of variant
        if len(ref) == len(cons):
            var_type = "SNP"
        elif len(ref) > len(cons):
            var_type = "DEL"
        else:
            var_type = "INS"
        # Add to variants
        msa_variants[anchor_pos] = [var_type, "II", anchor_nuc, ref, cons, freq, region, aa_mut]

        # Sort dictionary by position
        msa_variants_sorted = {key:msa_variants[key] for key in sorted(msa_variants)}

    with open(f"{output_dir}/{output_name}_anchored_variants.txt","w") as f:
        f.write("\t".join(["Anchor Position","Variant", "Variant Type", "Anchor Allele","Strain Allele", "Sample Allele",
            "Allele Frequency", "Protein","AA mutation"]))
        f.write("\n")
        for k in msa_variants_sorted.keys():
            f.write(str(k)+"\t"+"\t".join(msa_variants[k])+"\n")
    
