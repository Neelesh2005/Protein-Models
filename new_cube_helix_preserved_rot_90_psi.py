
import numpy as np
from Bio.PDB import PDBParser, PDBIO, PPBuilder
from Bio.PDB.vectors import Vector, rotaxis2m
from math import radians

def rotate_psi_preserve_helix(structure,chain_id, res_num, angle_deg):
    """
    Rotates the psi (ψ) angle of a given residue while preserving the helical structure.
    """

    
    target_res = None
    res_list = []
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                res_list = list(chain.get_residues())
                try:
                    res_idx = [r.id[1] for r in res_list].index(res_num)
                    target_res = res_list[res_idx]
                except (ValueError, KeyError):
                    raise ValueError(f"Residue {res_num} not found in chain {chain_id}")

    if not target_res:
        raise ValueError("Target residue not found")

    c_atom = target_res['C']
    next_res_index = res_list.index(target_res) + 1
    if next_res_index >= len(res_list):
        raise ValueError("Cannot rotate psi: No next residue available")
    
    next_res = res_list[next_res_index]
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
    return structure
    # io = PDBIO()
    # io.set_structure(structure)
    # output_file = f"helix_rotated_psi_{res_num}_{angle_deg}deg.pdb"
    # io.save(output_file)
    # print(f"Successfully rotated structure saved to {output_file}")
    
# def random_residue(pdb_file, chain_id):
#     """
#     Returns a random residue number from the given chain.
#     """
#     parser = PDBParser(QUIET=True)
#     structure = parser.get_structure("helix", pdb_file)
    
#     for model in structure:
#         for chain in model:
#             if chain.id == chain_id:
#                 residues = list(chain.get_residues())
#                 residue = residues[np.random.randint(1, len(residues))]  # Avoid first residue (no φ)
#                 return residue.id[1]

# try:
#     res_num = random_residue(r"e:\protein-model\8U1T_correct.pdb", "A")
#     rot_angle = np.random.randint(90,180)
#     rotate_psi_preserve_helix(
#         pdb_file=r"e:\protein-model\8U1T_correct.pdb",
#         chain_id="A",
#         res_num=res_num,
#         angle_deg=rot_angle
#     )
# except Exception as e:
#     print(f"Error: {str(e)}")
parser = PDBParser(QUIET=True)
structure = parser.get_structure("helix", r"e:\\protein-model\\8U1T_correct.pdb")
rotate_psi_preserve_helix(structure,"A", 14, 90)
