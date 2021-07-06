#!/usr/bin/env python3

import sys
import re

if len(sys.argv) != 4:
    print('ERROR: "Input format required is: <protein sequence> <original model file (xml format)> <new_model file (xml format)>"')
    sys.exit(1)

protein_seq = sys.argv[1]
original_model_file = sys.argv[2]
new_model_file = sys.argv[3]

codon_optimization ={
    'A': 'GCU', 'C': 'UGU', 'D': 'GAC', 'E': 'GAG', 'F': 'UUC',
    'G': 'GGU', 'H': 'CAC', 'I': 'AUU', 'K': 'AAG', 'L': 'UUG',
    'M': 'AUG', 'N': 'AAC', 'P': 'CCA', 'Q': 'CAA', 'R': 'AGA',
    'S': 'UCU', 'T': 'ACU', 'V': 'GUU', 'W': 'UGG', 'Y': 'UAC'
    }

RNA_seq = ''
for aa in protein_seq:
    RNA_seq += codon_optimization[aa]
#add the stop codon (UAA is the most likely for P.pastoris with 75%)
RNA_seq += 'UAA'


reverse_transcription = str.maketrans('AUGC', 'ATCG')
DNA_seq = RNA_seq.translate(reverse_transcription)

sbml_aa_to_id = {
    'A': 'M_m153', 'R': 'M_m158', 'N': 'M_m161', 'D': 'M_m156',
    'C': 'M_m331', 'E': 'M_m154', 'Q': 'M_m163', 'G': 'M_m210',
    'H': 'M_m490', 'I': 'M_m239', 'L': 'M_m227', 'K': 'M_m203',
    'M': 'M_m343', 'F': 'M_m272', 'P': 'M_m185', 'S': 'M_m279',
    'T': 'M_m305', 'W': 'M_m280', 'Y': 'M_m274', 'V': 'M_m222'
    }
    
sbml_id_to_aa = {id:aa for aa,id in sbml_aa_to_id.items()}

sbml_dna = {'M_m404': 'A', 'M_m437': 'T', 'M_m431': 'C', 'M_m389': 'G'}

sbml_rna = {'M_m94': 'A', 'M_m418': 'U', 'M_m423': 'C', 'M_m384': 'G'}

dna_MW = {'A': 313.2, 'T': 304.2, 'G': 329.2, 'C': 289.2}

rna_MW = {'A': 329.2, 'U': 306.2, 'G': 345.2, 'C': 305.2}

aa_MW = {
    'A': 71.078,  'R': 156.186, 'N': 114.103, 'D': 115.087, 'C': 103.143,
    'E': 128.129, 'Q': 129.114, 'G': 57.051,  'H': 137.139, 'I': 113.158,
    'L': 113.158, 'K': 128.172, 'M': 131.196, 'F': 147.174, 'P': 97.115,
    'S': 87.077,  'T': 101.104, 'W': 186.21,  'Y': 163.173, 'V': 99.131
    }

def get_stoichiometry(seq, sequence_type, MW_dict):
    seq_len = len(seq)
    counts = {}
    sequence_type = sequence_type.upper() # to allow user introduce lower and uppercase (should be fixed with app, because we control the input!)
    # Step 1: Get count for each amino acid 
    for bp in seq:
        if bp in counts:
            counts[bp] += 1
        else: 
            counts[bp] = 1
    # we need to consider the case that some protein might not contain some specific amino acids, and we need to specify that it is 0
    if sequence_type == 'PROTEIN':
        for amino_acid in sbml_aa_to_id.keys():
            if amino_acid not in counts.keys():
                counts[amino_acid] = 0

    # Step 2: Get gr/mol of protein for each amino acid
    gr_mol = {}
    for bp in counts:
        gr_mol[bp] = MW_dict[bp] * counts[bp] / seq_len

    # Step 4: Get mmol/ gr of protein for each amino acid
    mmol_gr = {}
    for bp in counts:
        mmol_gr[bp] = round(gr_mol[bp] / sum(gr_mol.values()) / MW_dict[bp] * 1000, 5)
    
    if sequence_type == 'DNA':
        atp_multiplier = 3.4
    elif sequence_type == 'RNA':
        atp_multiplier = 2.4
    elif sequence_type == 'PROTEIN':
        atp_multiplier = 4.3
    
    atp = sum(mmol_gr.values()) * atp_multiplier; mmol_gr['ATP'] = round(atp, 5)
    
    return (mmol_gr)

dna_stoch = get_stoichiometry(DNA_seq, 'dna', dna_MW)
rna_stoch = get_stoichiometry(RNA_seq, 'rna', rna_MW)
aa_stoch = get_stoichiometry(sys.argv[1], 'protein', aa_MW)
print(dna_stoch)
print(rna_stoch)
print(aa_stoch)


# and also add it to sbml_id_to_aa to ensure that when we read them on the file, we can translate the id to the proper species
sbml_dna['M_m1'] = 'ATP'; sbml_rna['M_m1'] = 'ATP'; sbml_id_to_aa['M_m1'] = 'ATP'

try:
    input_file = open(original_model_file, 'r')
    output_file = open(new_model_file, 'w')
except IOError as err:
    print("Cant open file:", str(err));
    sys.exit(1)

reaction_flag = False

#dict for cleaner coding
reaction_dict = {
    'R_r1100': [aa_stoch, sbml_id_to_aa],
    'R_r1101': [dna_stoch, sbml_dna],
    'R_r1103': [rna_stoch, sbml_rna]
    }

for line in input_file: 
    
    #copy file line by line while looking for reaction
    if not reaction_flag:
        output_file.write(line)
        
        #look for reaction ID's specified in dict
        for reaction in reaction_dict:
            if line.find('id="' + reaction) != -1:
                
                #when reaction lines found, update flag and dictionaries
                reaction_flag = True
                stochiometry_dict = reaction_dict[reaction][0]
                sbml_ids = reaction_dict[reaction][1]
    

    #when reaction lines are reached
    if reaction_flag:
        
        #regex group 1 is the species, group 2 is the stoichiometric value
                           #<speciesReference species="M_m1360" stoichiometry="0.997" constant="true"/>
        regex = re.search(r'<speciesReference species="(.*)" stoichiometry="(.*)" constant="true"/>', line)
        if regex:
            species = regex.group(1)
            old_stochiometry = regex.group(2)
            
            #M_m3, -5 and -7 are all equivalent to M_m1
            if species in ('M_m3', 'M_m5', 'M_m7'):
                species = 'M_m1'
            
            #don't edit lines with stoic. value of 1
            if old_stochiometry == "1":
                output_file.write(line)
            #rewrite line with new value
            else:
                new_stochiometry = str(stochiometry_dict[sbml_ids[species]])
                print(new_stochiometry)
                output_file.write(line.replace(old_stochiometry, new_stochiometry))
        
        #keep copying file if reaction lines not found
        else:
            output_file.write(line)
            
            #update flag if end of reaction lines
            if line.find('</listOfProducts>') != -1:
                reaction_flag = False

