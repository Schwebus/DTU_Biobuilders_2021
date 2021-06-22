import sys 

aa_sequence = sys.argv[1]
protein_len = len(aa_sequence)
proportions = {}

for aa in aa_sequence:
    if aa in proportions:
        proportions[aa] += 1
    else: 
        proportions[aa] = 1

for aa in proportions:
    proportions[aa] /= protein_len
    
print(proportions)
