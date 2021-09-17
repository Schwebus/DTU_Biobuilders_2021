#!/usr/bin/env python3

import sys

def transcription(sequence):
    DNA = 'ATCG'
    RNA = 'UAGC'
    transcription = sequence.maketrans(DNA, RNA)

    RNA_seq = sequence.translate(transcription)
    return RNA_seq

def get_stoichiometry(seq, sequence_type, protein_name):

    # define different species ids on the sbml model
    sbml_aa= {
        'A': 'ala__L_c', 'R': 'arg__L_c', 'N': 'asn__L_c', 'D': 'asp__L_c',
        'C': 'cys__L_c', 'E': 'glu__L_c', 'Q': 'gln__L_c', 'G': 'gly_c',
        'H': 'his__L_c', 'I': 'ile__L_c', 'L': 'leu__L_c', 'K': 'lys__L_c',
        'M': 'met__L_c', 'F': 'phe__L_c', 'P': 'pro__L_c', 'S': 'ser__L_c',
        'T': 'thr__L_c', 'W': 'trp__L_c', 'Y': 'tyr__L_c', 'V': 'val__L_c'
        }

    sbml_dna = {'A': 'damp_c', 'T': 'dtmp_c', 'C': 'dcmp_c', 'G': 'dgmp_c'}

    sbml_rna = {'A': 'amp_c', 'U': 'ump_c', 'C': 'cmp_c', 'G': 'gmp_c'}

    dna_MW = {'A': 331.2, 'T': 322.2, 'G': 347.2, 'C': 307.2}

    rna_MW = {'A': 347.2, 'U': 324.2, 'G': 363.2, 'C': 323.2}

    aa_MW = {
        'A': 89.1,  'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2,
        'E': 147.1, 'Q': 146.2, 'G': 75.1,  'H': 155.2, 'I': 131.2,
        'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 'P': 115.1,
        'S': 105.1,  'T': 119.1, 'W': 204.2,  'Y': 181.2, 'V': 117.1
        }

    seq_len = len(seq)
    counts = {}
    sequence_type = sequence_type.upper() # to allow user introduce lower and uppercase 
    # Step 1: Get count for each amino acid, or DNA/RNA base
    for bp in seq:
        if bp in counts:
            counts[bp] += 1
        else: 
            counts[bp] = 1
    # we need to consider the case that some protein might not contain some specific amino acids, and we need to specify that it is 0 and tell the MW dict to get the info from
    if sequence_type == 'AA':
        MW_dict = aa_MW
        model_ids = sbml_aa
        for amino_acid in sbml_aa.keys():
            if amino_acid not in counts.keys():
                counts[amino_acid] = 0
    elif sequence_type == 'DNA':
        MW_dict = dna_MW
        model_ids = sbml_dna
    elif sequence_type == 'RNA':
        MW_dict = rna_MW
        model_ids = sbml_rna
    else:
        print("ERROR: The sequence type introduced is not correct. Options are DNA, RNA or AA")
        sys.exit(1)
    
    # Step 2: Get gr/mol of protein for each amino acid
    gr_mol = {}
    for bp in counts:
        gr_mol[bp] = MW_dict[bp] * counts[bp] / seq_len

    # Step 3: Get mmol/ gr of protein for each amino acid
    mmol_gr = {}
    for bp in counts:
        mmol_gr[bp] = round(gr_mol[bp] / sum(gr_mol.values()) / MW_dict[bp] * 1000, 5)
    
    # change it to the ids used for the model, as we need it to automize the reaction definition
    mmol_gr_sbml = {}
    for bp in mmol_gr:   
        mmol_gr_sbml[model_ids[bp]] = - mmol_gr[bp]
    
    # we add the ATP stoichiometry, as well as ADP and water accordingly, involved in the DNA, RNA and protein synthesis
    if sequence_type == 'DNA':
        atp_multiplier = 3.4
    elif sequence_type == 'RNA':
        atp_multiplier = 2.4
    elif sequence_type == 'AA':
        atp_multiplier = 4.3

    atp = sum(mmol_gr.values()) * atp_multiplier; 

    mmol_gr['ATP'] = round(atp, 5); 
    for count, metabolite in enumerate(['atp_c', 'h2o_c', 'adp_c']):
        if count < 2:
            mmol_gr_sbml[metabolite] = - round(atp, 5)
        else:
            mmol_gr_sbml[metabolite] = round(atp, 5)
    
    # finally include the product
    mmol_gr_sbml[protein_name + '_' + sequence_type] = 1 

    return mmol_gr_sbml