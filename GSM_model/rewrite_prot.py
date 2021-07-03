import sys
import re
import fileinput
# MGAFTEKQEALVSSSFEAFKANIPQYSVVFYTSILEKAPAAKDLFSFLSNGVDPSNPKLTGHAEKLFGLVRDSAGQLKANGTVVADAALGSIHAQKAITDPQFVVVKEALLKTIKEAVGDKWSDELSSAWEVAYDELAAAIKKAF

aa_sequence = sys.argv[1]
model_file = sys.argv[2]

sbml_aa_to_id = {'A': 'M_m153' , 'R': 'M_m158', 'N': 'M_m161' , 'D': 'M_m156', 'C': 'M_m331', 'E': 'M_m154', 'Q': 'M_m163', 'G': 'M_m210',
'H': 'M_m490', 'I': 'M_m239', 'L': 'M_m227' , 'K': 'M_m203', 'M': 'M_m343', 'F': 'M_m272', 'P': 'M_m185', 'S': 'M_m279',
'T': 'M_m305', 'W': 'M_m280', 'Y': 'M_m274', 'V': 'M_m222'}

sbml_dna = {'M_m404': 'A', 'M_m437': 'T', 'M_m431': 'C', 'M_m389': 'G', 'M_m1': 'ATP', 'M_m5': 'H2O'}

sbml_rna = {'M_m94': 'A', 'M_m418': 'U', 'M_m423': 'C', 'M_m384': 'G', 'M_m1': 'ATP', 'M_m5': 'H2O'}

sbml_id_to_aa = {id:aa for aa,id in sbml_aa_to_id.items()}

# define the molecular weight for each nucleotide, ribonucleotide and aa (obtained from excel)
dna_MW = {'A': 313.2, 'T': 304.2, 'G': 329.2, 'C': 289.2}

rna_MW = {'A': 329.2, 'U': 306.2, 'G': 345.2, 'C': 305.2}

aa_MW = {'A': 71.078, 'R': 156.186, 'N': 114.103, 'D': 115.087, 'C': 103.143, 'E': 128.129, 'Q': 129.114, 'G': 57.051,
'H': 137.139, 'I': 113.158, 'L': 113.158, 'K': 128.172, 'M': 131.196, 'F': 147.174, 'P': 97.115, 'S': 87.077,
'T': 101.104, 'W': 186.21, 'Y': 163.173, 'V': 99.131}


protein_len = len(aa_sequence)
aa_counts = {}

# Step 1: Get count for each amino acid 
for aa in aa_sequence:
    if aa in aa_counts:
        aa_counts[aa] += 1
    else: 
        aa_counts[aa] = 1

# we need to consider the case that some protein might not contain some specific amino acids, and we need to specify that it is 0
for amino_acid in sbml_aa_to_id.keys():
    if amino_acid not in aa_counts.keys():
        aa_counts[amino_acid] = 0

# Step 2: Get gr/mol of protein for each amino acid
gr_mol = {}
for aa in aa_counts:
    gr_mol[aa] = aa_MW[aa] * aa_counts[aa] / protein_len

# Step 4: Get mmol/ gr of protein for each amino acid
mmol_gr = {}
for aa in aa_counts:
    mmol_gr[aa] = round(gr_mol[aa] / sum(gr_mol.values()) / aa_MW[aa] * 1000, 5)

atp = sum(mmol_gr.values()) * 4.3; mmol_gr['ATP'] = round(atp, 5)
# and also add it to sbml_id_to_aa to ensure that when we read them on the file, we can translate the id to the proper species
sbml_id_to_aa['M_m1'] = 'ATP'; 
print(mmol_gr)
reaction_flag = False
with fileinput.FileInput(model_file, inplace=True) as file:
    for line in file: 
        if line.find('##') != -1:
            reaction_flag = True 
        if not reaction_flag : 
            print(line, end = '')          
        if reaction_flag:
            if re.search('speciesReference', line):
                species = re.search(r'species="(.*)" stoichiometry', line).group(1) # get the id of the species for each reactant and product of the reaction
                old_stochiometry = re.search(r'stoichiometry="(.*)" constant', line).group(1)
                if old_stochiometry == "1":
                    new_stochiometry == "1"
                    print(line, end = '')
                elif species in ('M_m3', 'M_m5', 'M_m7'):
                    new_stochiometry = mmol_gr[sbml_id_to_aa['M_m1']]
                    print(line.replace(old_stochiometry, str(new_stochiometry)), end = '')
                else:
                    new_stochiometry = mmol_gr[sbml_id_to_aa[species]]
                    print(line.replace(old_stochiometry, str(new_stochiometry)), end = '')
            elif re.search('/listOfProducts',line):
                print(line, end = '')
                reaction_flag = False
            else:
                print(line, end = '')
 