from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.internal_coords import IC_Chain

# Load structure and generate internal coordinates
parser = PDBParser(QUIET=True)
structure = parser.get_structure("input", r"E:\protein-model\8U1T_modifi.pdb")

# Define target residue by (chain_id, res_id)
TARGET_CHAIN = "A"
TARGET_RESIDUE = (" ", 13, " ")  # Residue 5 in chain A

for model in structure:
    for chain in model:
        if chain.id == TARGET_CHAIN:
            ic_chain = IC_Chain(chain)
            ic_chain.atom_to_internal_coordinates(verbose=False)
            
            # Find target residue
            residues = list(chain.get_residues())
            try:
                target_index = residues.index(chain[TARGET_RESIDUE])
                if target_index == 0:
                    print("First residue has no phi angle")
                    continue
            except (KeyError, ValueError):
                print(f"Residue {TARGET_RESIDUE} not found in chain {TARGET_CHAIN}")
                continue

            # Modify phi angle
            target_res = chain[TARGET_RESIDUE]
            if target_res.internal_coord:
                current_phi = target_res.internal_coord.get_angle("phi")
                if current_phi is not None:
                    target_res.internal_coord.set_angle("phi", current_phi + 90)

            # Rebuild coordinates
            ic_chain.internal_to_atom_coordinates(verbose=False)

# Save modified structure
io = PDBIO()
io.set_structure(structure)
io.save("modified_phi_ILE.pdb")
