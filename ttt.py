from Bio import *
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
structure = PDBParser().get_structure('8U1T', r'E:\protein-model\8U1T_correct.pdb')    
import numpy as np

model = structure[0]
chain = model['A']

# for i in chain.get_residues():
#     print(type(i.get_resname()))

res_list = list(chain.get_residues())
length = len(res_list)
print(length)
print(res_list[0].get_id()[1]+length,res_list[0].get_id()[1])
low = res_list[0].get_id()[1]
high = res_list[0].get_id()[1]+length
print(np.random.randint(low,high))
# sequence_name = res_list[np.random.randint(res_list[0].get_id()[1],res_list[0].get_id()[1]+length)].get_resname()
# print(sequence_name)