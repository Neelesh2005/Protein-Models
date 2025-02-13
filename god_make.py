import numpy as np
from Bio.PDB import PDBParser, PDBIO, DSSP
from Bio.PDB.vectors import Vector, rotaxis2m
from math import radians
import os

def rotate_phi_with_helix(pdb_file, chain_id, res_num, angle_deg):
    """
    Rotate phi angle while preserving helical structure and secondary annotations
    """
    # Parse structure with preserved header information
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("rotated", pdb_file)
    
    # Preserve original header and secondary structure records
    with open(pdb_file) as f:
        header = [line for line in f if line.startswith(("HELIX","SHEET","TURN"))]

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
                    
                    # Collect all atoms from CA onward (including downstream residues)
                    downstream_atoms = []
                    for res in residues[res_idx:]:
                        for atom in res:
                            if res == target_res and atom.name in ['N', 'CA']:
                                continue  # Keep N/CA fixed as rotation pivot
                            downstream_atoms.append(atom)
                except ValueError:
                    raise ValueError(f"Residue {res_num} not found in chain {chain_id}")

    if not target_res:
        raise ValueError("Target residue not found")

    # Get rotation axis (N-CA vector)
    n_coord = target_res['N'].coord
    ca_coord = target_res['CA'].coord
    axis = Vector(ca_coord - n_coord).normalized()

    # Create rotation matrix using BioPython's optimized function
    angle_rad = radians(angle_deg)
    rotation_matrix = rotaxis2m(angle_rad, axis)

    # Apply rotation to downstream atoms
    for atom in downstream_atoms:
        translated = atom.coord - n_coord
        rotated = np.dot(translated, rotation_matrix.T)
        atom.coord = rotated + n_coord

    # Save with preserved secondary structure records
    output_file = f"pleasdeg.pdb"
    
    # Write header first
    with open(output_file, "w") as f:
        f.writelines(header)
        
        # Then write atomic coordinates
        io = PDBIO()
        io.set_structure(structure)
        io.save(f, write_end=False)  # Prevent END record
        
    print(f"Saved rotated structure to {output_file}")
    
    # Recalculate secondary structure with DSSP
    print("\nRecomputing secondary structure with DSSP:")
    model = structure[0]
    dssp = DSSP(
    model, 
    output_file,
    dssp=r'C:\Miniconda3\Library\bin\mkdssp.exe'  # Full path to executable
)

    
    # Update HELIX records based on DSSP results
    helix_residues = []
    for key in dssp.keys():
        if dssp[key][2] == 'H':  # Helix residues
            helix_residues.append(key[0][1])
    
    print(f"Helix residues after rotation: {helix_residues}")

# Example usage
if __name__ == "__main__":
    input_pdb = r"E:\protein-model\8U1T_correct.pdb" 
    rotate_phi_with_helix(
        pdb_file=input_pdb,
        chain_id="A",
        res_num=33,  # ILE A13 from sample data
        angle_deg=90
    )
