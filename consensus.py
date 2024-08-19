import os
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
def iupac_code(nts):
    iupac_code_dict = {
        'AC':'M',
        'AG':'R',
        'AT':'W',
        'CG':'S',
        'CT':'Y',
        'GT':'K',
        'ACG':'V',
        'ACT':'H',
        'AGT':'D',
        'CGT':'B',
        'ACGT':'N'
    }
    rev_iupac_code = {iupac_code_dict[k]:k for k in iupac_code_dict.keys()}
    corrected_nt = []

    for nt in nts:
        if nt not in ["A", "C", "G", "T"]:
            corrected_nt += [e for e in rev_iupac_code[nt]]  ## if nuc is amb, then change to double nucleotide (CG)
        else:
            corrected_nt += [nt]

    unique_nts = "".join((sorted(list(np.unique(corrected_nt)))))

    return iupac_code_dict[unique_nts]


def extract_headers_and_sequences(ref):
    header = ''
    sequences = {}
    for line in open(ref,'r'):
        if '>' in line:
            # initialize dict entry sequences[header] = sequence
            header = line.replace('\n','').replace('>','').split(' ')[0]
            sequences[header] = ''
        else:
            sequences[header] += line.replace('\n','')
    return sequences

def apply_variants(sequences,vcf):
    consensus_sequences = dict(sequences)
    corrected_position_index = 0
    prev_chrom = ""
    vcf_lines = []
    for line in open(vcf,'r'):
        if '#' in line:
            pass
        else:
            vcf_lines.append(line)
    i = 0
    while i < len(vcf_lines):
        indel = False
        line = vcf_lines[i].replace('\n','').split('\t')
        chrom = line[0]
        pos = int(line[1])
        ref = [line[3]]
        alt = [line[4]]
        # print(chrom,pos,ref,alt)
        # Check if new chromosome, if so, reset corrected position index
        if chrom != prev_chrom:
            corrected_position_index = 0
        if len([e for e in ref if len(e)>1]) > 0 or len([e for e in alt if len(e)>1]) > 0:
            indel = True
        af= [float(line[7].split(';')[1].replace('AF=',''))]
        multialleles=True
        ii = 1
        while multialleles==True:
            if i+ii < len(vcf_lines):
                next_line = vcf_lines[i+ii].replace('\n','').split('\t')
                next_pos = int(next_line[1])
                if next_pos == pos:
                    if len(next_line[4]) - len(next_line[3]) != 0:
                        indel = True
                    ref.append(next_line[3])
                    alt.append(next_line[4])
                    af.append(float(next_line[7].split(';')[1].replace('AF=','')))
                    ii+=1
                else:
                    multialleles=False
            else:
                multialleles=False
        # print(ref,alt,af,i,ii)
        # Take max allele frequency for each indel to determine
        # if variant is a minority indel, pass
        if indel == True and max([af[j] for j in range(len(af)) if len(alt[j])-len(ref[j]) != 0]) < 0.5:
            pass
        # If there is only one alternate allele
        elif len(af) == 1:
            ref = ref[0]
            alt = alt[0]
            af = af[0]
            # Calculate variant length
            variant_length = len(alt) - len(ref)
            # If the variant is a SNP
            if variant_length == 0 and len(ref) == 1:
                # Calculate correct positions to split string
                before = pos-1+corrected_position_index
                after = pos+corrected_position_index
                # If the variant is a major SNP, change ref nt to alt nt
                if af > 0.95:
                    consensus_sequences[chrom] = consensus_sequences[chrom][:before] + alt + consensus_sequences[chrom][after:]
                # If the variant is a minority SNP, then find appropriate iupac ambiguity code and insert
                else:
                    amb = iupac_code(''.join(sorted([ref,alt])))
                    consensus_sequences[chrom] = consensus_sequences[chrom][:before] + amb + consensus_sequences[chrom][after:]
            # If the variant is an indel
            else:
                #it's AF is >= 0.5 because of previous filter
                if af > 0.5:
                    # calculate corrected positions to split string
                    before = pos-1+corrected_position_index
                    after = pos-1+len(ref)+corrected_position_index
                    consensus_sequences[chrom] = consensus_sequences[chrom][:before] + alt + consensus_sequences[chrom][after:]
                    corrected_position_index += len(alt)-len(ref)
                else:
                    pass
        # If there are multiple alternate alleles
        else:
            # If all alternate alleles are SNPs
            if len([e for e in alt if len(e)>1]) == 0 and len([e for e in ref if len(e)>1]) == 0:
                before = pos-1+corrected_position_index
                after = pos+corrected_position_index
                amb = iupac_code(''.join(sorted(list(set(ref+alt)))))
                consensus_sequences[chrom] = consensus_sequences[chrom][:before] + amb + consensus_sequences[chrom][after:]
            # Apply allele with greatest frequency
            else:
                # Caclulate max af and associated alternate allele
                # It must be <= 0.95 because of previous check
                max_af = 0
                max_index = 0
                for j in range(len(af)):
                    if af[j] > max_af:
                        max_af = af[j]
                        max_index = j
                alt = alt[max_index]
                ref = ref[max_index]
                # Calculate variant length
                variant_length = len(alt) - len(ref)
                # If the variant is a SNP
                if variant_length == 0 and len(ref) == 1:
                    # calculate corrected positions to split string
                    before = pos-1+corrected_position_index
                    after = pos+corrected_position_index
                    # The variant AF can't be over 0.95, therefore find appropriate iupac ambiguity code and insert
                    amb = iupac_code(''.join(sorted([ref,alt])))
                    consensus_sequences[chrom] = consensus_sequences[chrom][:before] + amb + consensus_sequences[chrom][after:]
                # If the variant is a substitution, it's af is >= 0.5 because of previous filtering
                elif variant_length == 0 and len(ref) > 1:
                    # calculate corrected positions to split string
                    before = pos-1+corrected_position_index
                    after = pos+len(ref)+corrected_position_index
                    consensus_sequences[chrom] = consensus_sequences[chrom][:before] + alt + consensus_sequences[chrom][after:]
                # If the variant is an insertion, it's af is >= 0.5 because of previous filtering
                elif len(alt) > len(ref):
                    # calculate corrected positions to split string
                    before = pos-1+corrected_position_index
                    after = pos+corrected_position_index
                    consensus_sequences[chrom] = consensus_sequences[chrom][:before] + alt + consensus_sequences[chrom][after:]
                    corrected_position_index += len(alt)-len(ref)
                # If the variant is a deletion, it's af is >= 0.5 because of previous filtering
                elif len(ref) > len(alt):
                    # calculate corrected positions to split string
                    before = pos+corrected_position_index
                    after = pos+len(ref)-len(alt)+corrected_position_index
                    consensus_sequences[chrom] = consensus_sequences[chrom][:before] + consensus_sequences[chrom][after:]
                    corrected_position_index += len(alt)-len(ref)
        i += ii
        prev_chrom = chrom
    return consensus_sequences

def mask(sequences,bed):
    masked_consensus_sequences = dict(sequences)
    for line in open(bed,'r'):
        line = line.replace('\n','').split('\t')
        chrom = line[0]
        before = int(line[1])
        after = int(line[2])
        length = after-before
        masked_consensus_sequences[chrom] = masked_consensus_sequences[chrom][:before] + 'N'*length + masked_consensus_sequences[chrom][after:]
    return masked_consensus_sequences


def main():
    parser = argparse.ArgumentParser(description='Generate a consensus sequence that contains ambiguous variants. \
    \n usage: python consensus.py -r reference.fasta -v variants.vcf -b low_cov.bed -o output.fasta',
    formatter_class=RawTextHelpFormatter)
    parser.add_argument('-r', metavar='reference',
     type=str, nargs='?', required=True, help='The path to the reference fasta file')
    parser.add_argument('-v', metavar='variants',
     type=str, nargs='?', required=True, help='The path to the variants vcf file')
    parser.add_argument('-b', metavar='low coverage regions',
     type=str, nargs='?', required=False, help='The path to the low coverage regions bed file')
    parser.add_argument('-o', metavar='output',
     type=str, nargs='?', required=True, help='The desired output filename')
    args = parser.parse_args()
    print("Generating consensus.")
    if os.path.isfile(args.r):
        if os.path.isfile(args.v):
            print("Pulling out reference sequence")
            sequences = extract_headers_and_sequences(args.r)
            print("Applying variants")
            consensus_sequences = apply_variants(sequences, args.v)
            # bed file provided
            if args.b:
                if os.path.isfile(args.b):
                    print("Masking sequence")
                    masked_consensus_sequences = mask(consensus_sequences, args.b)
                    with open(args.o,'w') as f:
                        for key in masked_consensus_sequences.keys():
                            f.write('>'+'consensus_'+key)
                            f.write('\n')
                            f.write(masked_consensus_sequences[key])
                            f.write('\n')
                else:
                    print('Invalid bed file path.')
            # bed file not provided
            else:
                with open(args.o,'w') as f:
                    for key in consensus_sequences.keys():
                        f.write('>'+'consensus_'+key)
                        f.write('\n')
                        f.write(consensus_sequences[key])
                        f.write('\n')
        else:
            print('Invalid vcf file path.')
    else:
        print("Invalid reference fasta path.")


if __name__ == "__main__":
    main()
