## INPUTS: GFF, reference .fasta, and VCF
## Output: amino acid mutations
## Optional output: visualization of mutations
import os
import sys
# from BCBio import GFF
from dna_features_viewer import GraphicFeature, GraphicRecord
from dna_features_viewer import BiopythonTranslator
from Bio import SeqIO
from bokeh.resources import CDN
from bokeh.embed import file_html
# from Bio import pairwise2


print("Command:",' '.join(sys.argv))

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

revcom = {"A":"T","T":"A","G":"C","C":"G"}

variants = {}
for line in open(sys.argv[1]):
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

if len(list(variants.keys())) == 0:
    print("No Variants in VCF")
    quit()


# for line in open(sys.argv[2]):
#     if '>' not in line:
#         reference += line.strip()
reference = ""
graphic_features = []
features = {}
for r in SeqIO.parse(open(sys.argv[3],"r"), "genbank"):
    if reference == "":
        reference = ''.join(list(r.seq))
    for f in r.features:
        # print("Type: ", f.type)
        if f.type == "source":
            org_id = sys.argv[3].replace(".gb","")
            if 'organism' in f.qualifiers:
                org_id += f" {f.qualifiers['organism'][0]}"
            if 'strain' in f.qualifiers:
                org_id += f" {f.qualifiers['strain'][0]}"
            source = GraphicFeature(start=f.location.start,end=f.location.end,strand=+1,label=org_id,data=org_id,color="#f9cb9c")
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
                graphic_features.append(GraphicFeature(start=f.location.start,end=f.location.end,strand=f.location.strand,label=label,data=label,color="#8e7cc3"))
            ### For ensuring correct mucleotide sequence and translation
            # aas = ''.join([genetic_code(seq[i:i+3].upper()) for i in range(0,len(seq),3)])
            # if aas[:-1] == f.qualifiers['translation'][0]:
                # print("correct translation")

def get_mutated_region(region,reference,position,ref,alt):
    # set for 0 index used in biopython
    position -= 1
    # If deletion, then store deleted positions
    deletion_positions = []
    if len(ref) > len(alt):
        for i in range(1,len(ref)):
            deletion_positions.append(position+i)
    # Read genbank and pull out sequence of region adding in the variant
    for r in SeqIO.parse(open(sys.argv[3],"r"), "genbank"):
        for f in r.features:
            if f.type == 'CDS':
                try:
                    label = f.qualifiers['product']
                except:
                    label = f.qualifiers['note']
                # print("Supposed product: ", f.qualifiers['product'])
#            if f.type == 'CDS' and f.qualifiers['product'][0] == region:
                # print(label, region)
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

# Go through variants
for position in variants.keys():
    full_position = position
    # "_" is present if there are multiple variants at a position
    if "_" in position:
        position = position[:position.index("_")]
    # Find if variant is in CDS region
    # if int(position) == 13386 and len([e for e in variants[position][1][1:] if e != 'C']) == 0:
    #     variants[full_position].append("ORF1ab")
    #     variants[full_position].append("-1 Ribosomal frameshift site")
    if int(position) == -10:
        print("This should never happen. It's just a chunk of code to keep things working, since the above was commented out.")
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
                # if features[k]['strand'] == -1:
                #     region = ''.join([revcom[e] for e in nt_seq])
                # If the lengths of the ref and alt alleles are the same and == 1 then SNP
                if len(variants[position][0])==len(variants[position][1]) and len(variants[position][1]) == 1:
                    # print(protein,reference,int(position),variants[position][0],variants[position][1])
                    mutated_region = get_mutated_region(protein,reference,int(position),variants[position][0],variants[position][1])
                    if mutated_region == region:
                        print("mutated and region are same")
                        variants[full_position].append('UTR')
                        variants[full_position].append("n/a")
                    else:
                        # find the variant position relative to the region
                        # region_variant_position = int(position) - 1 - start
                        # if features[k]['strand'] == -1:
                        #     region_variant_position = int(position) - 1 - end    
                        # determine the codon number and codon position [0,1,2]
                        # codon_number = region_variant_position//3
                        # codon_position = region_variant_position - (codon_number*3) # 0, 1, or 2
                        # pull out the original and mutated codon and translate
                        # original_codon = region[(codon_number*3):(codon_number*3)+3]
                        # mutated_codon = mutated_region[(codon_number*3):(codon_number*3)+3]
                        # print(variants[position])
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
                    mutated_region = get_mutated_region(protein,reference,int(position),variants[position][0],variants[position][1])
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
                    mutated_region = get_mutated_region(protein,reference,int(position),variants[position][0],variants[position][1])
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
variant_file = sys.argv[1].replace(".vcf","")
with open(f"{variant_file}_variants.txt","w") as f:
    f.write("\t".join(["Position","Reference allele","Alternate allele","Allele frequency","Protein","AA mutation"]))
    f.write("\n")
    for k in variants.keys():
        f.write(k+"\t"+"\t".join(variants[k])+"\n")


######## Write to low_coverage.txt ##########
lowcov_file = sys.argv[1].replace("lofreq.vcf","")
with open(f"{variant_file}_lowcoverage.txt","w") as f:
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



###################
##################### Comment the plotting out for now for Lassa###################
#################
########## Plot variants ############
# Add genbank features and variant features to a list of features
feature_list = [source]
for gf in graphic_features:
    feature_list.append(gf)
for k in variants.keys():
    if "_" in k:
        pos = k[:-2]
    else:
        pos = k
    if variants[k][4] != "n/a":
        label = variants[k][4]
        html = variants[k][2]
    else:
        label = variants[k][0] + pos + variants[k][1]
        html = "intergenic "+ variants[k][2]
    if "aa" in label:
        label = "*"
    feature_list.append(GraphicFeature(start=int(pos),end=int(pos)+1,strand=+1,label=label,data=label,color="#e06666",html=html,box_color='#e06666',label_link_color='#e06666'))
# create graphic record for plotting and generate html
record = GraphicRecord(sequence=reference,features=feature_list)
plot = record.plot_with_bokeh(figure_width=8)
# create base html page 
output_name = sys.argv[1].replace(".vcf","")
with open(f"{output_name}_raw.html", "w+") as f:
    f.write(file_html(plot, CDN, "Example Sequence"))
# edit base html page to have title and figure description
with open(f"{output_name}.html", "w+") as f:
    for line in open(f"{output_name}_raw.html", "r"):
        f.write(line)
        if "<body>" in line:
            base_output = output_name.split("/")[0]
            f.write(f"<h3>Interactive display of {base_output} mutations relative to {org_id} annotations.<h3>\n")
        elif "</div>" in line:
            f.write("<p style=\"font-size:12px;font-weight:normal\">* indicates a significant protein mutation. See variants.txt for mutation details<p>\n")
# remove base html page
os.system(f"rm {output_name}_raw.html")
print("Completed") 
