from Bio.PDB import *
from Bio.PDB.vectors import Vector, rotaxis2m
import numpy as np
from math import radians

def rotate_residue(pdb_file, chain_id, res_num, phi_angle):
    """
    Rotate a single residue while maintaining structure connectivity
    """
    # Parse structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    
    # Get the specific residue
    model = structure[0]
    chain = model[chain_id]
    residue = chain[(' ', res_num, ' ')]
    
    # Get rotation axis (N-CA bond)
    n_atom = residue['N']
    ca_atom = residue['CA']
    
    # Calculate rotation axis as a Vector
    n_coord = n_atom.get_coord()
    ca_coord = ca_atom.get_coord()
    rotation_axis = Vector(ca_coord - n_coord)
    rotation_axis.normalize()
    
    # Create rotation matrix using Vector object
    rotation = rotaxis2m(radians(phi_angle), rotation_axis)
    
    # Apply rotation only to atoms after CA
    for atom in residue:
        if atom.name not in ['N', 'CA']:
            atom.transform(rotation, np.array([0., 0., 0.]))
    
    # Save modified structure
    io = PDBIO()
    io.set_structure(structure)
    output_file = f"rotated_{res_num}_{phi_angle}.pdb"
    io.save(output_file)
    print(f"Saved modified structure to {output_file}")

# Example usage
if __name__ == "__main__":
    input_pdb = r"E:\protein-model\8U1T_coords.pdb"
    chain_id = "A"
    residue_number = 13
    rotation_angle = 90
    
    rotate_residue(input_pdb, chain_id, residue_number, rotation_angle)
