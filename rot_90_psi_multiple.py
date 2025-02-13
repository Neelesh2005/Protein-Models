import numpy as np
from Bio.PDB import PDBParser, PDBIO, PPBuilder
from Bio.PDB.vectors import Vector, rotaxis2m
from math import radians

def rotate_psi_preserve_helix(pdb_file, chain_id, res_nums, angle_deg=90):
    """
    Rotates the psi (ψ) angle of given residues while preserving the helical structure.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("helix", pdb_file)
    
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                res_list = list(chain.get_residues())
                res_indices = {r.id[1]: idx for idx, r in enumerate(res_list)}
                
                for res_num in res_nums:
                    if res_num not in res_indices:
                        print(f"Residue {res_num} not found in chain {chain_id}")
                        continue
                    
                    target_res = res_list[res_indices[res_num]]
                    next_res_index = res_indices[res_num] + 1
                    if next_res_index >= len(res_list):
                        print(f"Cannot rotate psi for residue {res_num}: No next residue available")
                        continue
                    
                    next_res = res_list[next_res_index]
                    c_atom = target_res['C']
                    n_next_atom = next_res['N']
                    rotation_center = c_atom.coord
                    rotation_axis = Vector(n_next_atom.coord - c_atom.coord).normalized()
                    rotation_matrix = rotaxis2m(radians(angle_deg), rotation_axis)
                    
                    downstream_atoms = []
                    for residue in res_list[next_res_index:]:
                        for atom in residue:
                            downstream_atoms.append(atom)
                    
                    for atom in downstream_atoms:
                        translated = atom.coord - rotation_center
                        rotated = np.dot(rotation_matrix, translated)
                        atom.coord = rotated + rotation_center
                    
                    print(f"Rotated psi angle of residue {res_num} by {angle_deg} degrees")
    
    io = PDBIO()
    io.set_structure(structure)
    output_file = f"helix_rotated_psi_{'_'.join(map(str, res_nums))}_{angle_deg}deg.pdb"
    io.save(output_file)
    print(f"Successfully rotated structure saved to {output_file}")

def random_residues(pdb_file, chain_id, num_residues=3):
    """
    Returns a list of random residue numbers from the given chain.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("helix", pdb_file)
    
    res_nums = []
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                res_list = list(chain.get_residues())
                res_indices = list(range(len(res_list)))
                res_indices.remove(0)
                
                for _ in range(num_residues):
                    res_idx = np.random.choice(res_indices)
                    res_nums.append(res_list[res_idx].id[1])
                    res_indices.remove(res_idx)
    
    return res_nums
try:
    len_res = np.random.randint(1, 8)
    res_nums = random_residues(r"e:\\protein-model\\8U1T_correct.pdb","A",len_res)  # Example residue numbers
    rotate_psi_preserve_helix(
        pdb_file=r"e:\\protein-model\\8U1T_correct.pdb",
        chain_id="A",
        res_nums=res_nums,
        angle_deg=90
    )
except Exception as e:
    print(f"Error: {str(e)}")
