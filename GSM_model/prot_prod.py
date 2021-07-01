import sys 

# define the molecular weight for each nucleotide, ribonucleotide and aa (obtained from excel)
dna_MW = {'A': 313.2, 'T': 304.2, 'G': 329.2, 'C': 289.2}

rna_MW = {'A': 329.2, 'U': 306.2, 'G': 345.2, 'C': 305.2}

aa_MW = {'A': 71.078, 'R': 156.186, 'N': 114.103, 'D': 115.087, 'C': 103.143, 'E': 128.129, 'Q': 129.114, 'G': 57.051,
'H': 137.139, 'I': 113.158, 'L': 113.158, 'K': 128.172, 'M': 131.196, 'F': 147.174, 'P': 97.115, 'S': 87.077,
'T': 101.104, 'W': 186.21, 'Y': 163.173, 'V': 99.131}

aa_sequence = sys.argv[1]
protein_len = len(aa_sequence)
aa_counts = {}

# Step 1: Get count for each amino acid 
for aa in aa_sequence:
    if aa in aa_counts:
        aa_counts[aa] += 1
    else: 
        aa_counts[aa] = 1

# Step 2 (could be simplified with 3): Get amino acid percentage presence in the protein
aa_perc = {}
for aa in aa_counts:
    aa_perc[aa] = aa_counts[aa] / protein_len * 100

# Step 3: Get gr/mol of protein for each amino acid
gr_mol = {}
for aa in aa_counts:
    gr_mol[aa] = aa_MW[aa] * aa_perc[aa] / 100

# Step 4: Get mmol/ gr of protein for each amino acid
mmol_gr = {}
for aa in aa_counts:
    mmol_gr[aa] = round(gr_mol[aa] / sum(gr_mol.values()) / aa_MW[aa] * 1000, 5)

atp = sum(mmol_gr.values()) * 4.3

sbml_aa_nom = {'A': 153 , 'R': 158, 'N': 161 , 'D': 156, 'C': 331, 'E': 154, 'Q': 163, 'G': 210,
'H': 490, 'I': 239, 'L': 227 , 'K': 203, 'M': 343, 'F': 272, 'P': 185, 'S': 279,
'T': 305, 'W': 280, 'Y': 274, 'V': 222}

print(sbml_aa_nom)
print(mmol_gr)