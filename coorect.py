import numpy as np
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.vectors import Vector, rotaxis2m

def rotate_phi_with_downstream(pdb_file, chain_id, res_num, angle_deg):
    """
    Rotate phi angle of target residue and maintain connectivity
    by rotating downstream residues as rigid body
    
    Args:
        pdb_file: Input PDB file path
        chain_id: Chain identifier (str)
        res_num: Target residue number (int)
        angle_deg: Rotation angle in degrees (float)
    """
    
    # Parse structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("rotated", pdb_file)
    
    # Find target residue and downstream atoms
    target_res = None
    downstream_atoms = []
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                residues = list(chain.get_residues())
                try:
                    res_idx = [r.id[1] for r in residues].index(res_num)
                    target_res = residues[res_idx]
                    
                    # Collect downstream atoms (from CA onward)
                    downstream_atoms = []
                    for res in residues[res_idx:]:
                        for atom in res:
                            if res == target_res and atom.name in ['N', 'CA']:
                                continue  # Skip N and CA of target
                            downstream_atoms.append(atom)
                except ValueError:
                    raise ValueError(f"Residue {res_num} not found in chain {chain_id}")

    if not target_res:
        raise ValueError("Target residue not found")

    # Get rotation axis (N-CA vector)
    n_coord = target_res['N'].coord
    ca_coord = target_res['CA'].coord
    axis = Vector(ca_coord - n_coord).normalized()

    # Create rotation matrix
    angle_rad = np.radians(angle_deg)
    rotation_matrix = rotaxis2m(angle_rad, axis)

    # Apply rotation to downstream atoms
    for atom in downstream_atoms:
        translated = atom.coord - n_coord
        rotated = np.dot(translated, rotation_matrix.T)
        atom.coord = rotated + n_coord

    # Save modified structure
    io = PDBIO()
    io.set_structure(structure)
    output_file = f"rotated_phi_{res_num}_{angle_deg}deg.pdb"
    io.save(output_file)
    print(f"Saved rotated structure to {output_file}")

# Example usage:
if __name__ == "__main__":
    rotate_phi_with_downstream(
        pdb_file=r"E:\protein-model\8U1T_modifi.pdb",
        chain_id="A",
        res_num=33,  # ILE A13 from your sample data
        angle_deg=90
    )
