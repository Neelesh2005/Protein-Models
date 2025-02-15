import numpy as np
import logging
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.vectors import Vector, rotaxis2m
from math import radians

# Setup Logging
logging.basicConfig(level=logging.INFO, format="%(message)s")

def rotate_phi_preserve_helix(pdb_file, chain_id, res_num, angle_deg):
    """
    Rotates the phi (φ) angle of a given residue while preserving the helical structure.
    """
    logging.info(f"Loading PDB structure from: {pdb_file}")
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("helix", pdb_file)
    
    target_res = None
    chain = None

    for model in structure:
        for ch in model:
            if ch.id == chain_id:
                chain = ch
                residues = list(chain.get_residues())
                try:
                    res_idx = [r.id[1] for r in residues].index(res_num)
                    target_res = residues[res_idx]
                except ValueError:
                    raise ValueError(f"Residue {res_num} not found in chain {chain_id}")

    if not target_res:
        raise ValueError("Target residue not found")

    logging.info(f"Target residue: {target_res.resname} {res_num}")

    res_list = list(chain.get_residues())
    target_index = res_list.index(target_res)
    if target_index == 0:
        raise ValueError("Cannot rotate phi: No previous residue available")

    prev_res = res_list[target_index - 1]

    try:
        n_atom = target_res['N']
        ca_atom = target_res['CA']
    except KeyError:
        raise ValueError(f"Residue {res_num} does not have both N and CA atoms (possible issue with missing atoms).")

    rotation_center = n_atom.coord  # Rotate around N
    rotation_axis = Vector(ca_atom.coord - n_atom.coord).normalized()
    rotation_matrix = rotaxis2m(radians(angle_deg), rotation_axis)

    downstream_atoms = []
    for residue in res_list[target_index:]:
        for atom in residue:
            downstream_atoms.append(atom)

    logging.info(f"Rotating {len(downstream_atoms)} atoms about φ angle by {angle_deg} degrees")

    for atom in downstream_atoms:
        translated = atom.coord - rotation_center
        rotated = np.dot(rotation_matrix, translated)
        atom.coord = rotated + rotation_center

    io = PDBIO()
    io.set_structure(structure)
    output_file = f"helix_rotated_phi_{res_num}_{angle_deg}deg.pdb"
    io.save(output_file)
    logging.info(f"Saved rotated structure to {output_file}")

def random_residue(pdb_file, chain_id):
    """
    Returns a random residue number from the given chain.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("helix", pdb_file)
    
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                residues = list(chain.get_residues())
                residue = residues[np.random.randint(1, len(residues))]  # Avoid first residue (no φ)
                return residue.id[1]

# Example usage
try:
    res_num = random_residue(r"e:\protein-model\8U1T_correct.pdb", "A")
    rot_angle = np.random.randint(90,180)

    rotate_phi_preserve_helix(
        pdb_file=r"e:\protein-model\8U1T_correct.pdb",
        chain_id="A",
        res_num=res_num,
        angle_deg=rot_angle
    )
except Exception as e:
    logging.error(f"Error: {str(e)}")
