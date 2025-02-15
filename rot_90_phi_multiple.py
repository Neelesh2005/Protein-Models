import numpy as np
import logging
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.vectors import Vector, rotaxis2m
from math import radians

# Setup Logging
logging.basicConfig(level=logging.INFO, format="%(message)s")

def rotate_phi_preserve_helix(pdb_file, chain_id, res_nums, angle_deg=90):
    """
    Rotates the phi (φ) angle of given residues while preserving the helical structure.
    """
    logging.info(f"Loading PDB structure from: {pdb_file}")
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("helix", pdb_file)
    
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                res_list = list(chain.get_residues())
                res_indices = {r.id[1]: idx for idx, r in enumerate(res_list)}
                
                for res_num in res_nums:
                    if res_num not in res_indices:
                        logging.warning(f"Residue {res_num} not found in chain {chain_id}")
                        continue
                    
                    target_res = res_list[res_indices[res_num]]
                    target_index = res_indices[res_num]
                    
                    if target_index == 0:
                        logging.warning(f"Cannot rotate phi for residue {res_num}: No previous residue available")
                        continue
                    
                    prev_res = res_list[target_index - 1]
                    
                    try:
                        n_atom = target_res['N']
                        ca_atom = target_res['CA']
                    except KeyError:
                        logging.warning(f"Residue {res_num} does not have both N and CA atoms (possible issue with missing atoms).")
                        continue
                    
                    rotation_center = n_atom.coord
                    rotation_axis = Vector(ca_atom.coord - n_atom.coord).normalized()
                    rotation_matrix = rotaxis2m(radians(angle_deg), rotation_axis)
                    
                    # Collect all downstream atoms
                    downstream_atoms = []
                    for residue in res_list[target_index:]:
                        for atom in residue:
                            downstream_atoms.append(atom)
                    
                    logging.info(f"Rotating {len(downstream_atoms)} atoms about φ angle of residue {res_num} by {angle_deg} degrees")
                    
                    for atom in downstream_atoms:
                        translated = atom.coord - rotation_center
                        rotated = np.dot(rotation_matrix, translated)
                        atom.coord = rotated + rotation_center
    
    # Save modified structure
    io = PDBIO()
    io.set_structure(structure)
    output_file = f"helix_rotated_phi_{'_'.join(map(str, res_nums))}_{angle_deg}deg.pdb"
    io.save(output_file)
    logging.info(f"Saved rotated structure to {output_file}")


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
                res_indices.remove(0)  # Exclude first residue
                
                for _ in range(num_residues):
                    res_idx = np.random.choice(res_indices)
                    res_nums.append(res_list[res_idx].id[1])
                    res_indices.remove(res_idx)
    
    return res_nums
try:
    len_res = np.random.randint(1, 6)
    res_nums = random_residues(r"e:\\protein-model\\8U1T_correct.pdb","A",len_res)
    rot_angle = np.random.randint(90,180)
# Example residue numbers  # Example residue numbers
    rotate_phi_preserve_helix(
        pdb_file=r"e:\\protein-model\\8U1T_correct.pdb",
        chain_id="A",
        res_nums=res_nums,
        angle_deg=rot_angle
    )
except Exception as e:
    logging.error(f"Error: {str(e)}")
