from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.internal_coords import IC_Chain

parser = PDBParser(QUIET=True)
structure = parser.get_structure("input", r"E:\protein-model\8U1T_modifi.pdb")  # Raw string

for model in structure:
    for chain in model:
        ic_chain = IC_Chain(chain)
        ic_chain.atom_to_internal_coordinates(verbose=False)
        
        for residue in list(chain.get_residues())[1:]:  # Skip first residue
            if residue.internal_coord:
                current_phi = residue.internal_coord.get_angle("phi")
                if current_phi is not None:
                    residue.internal_coord.set_angle("phi", current_phi + 90)
        
        ic_chain.internal_to_atom_coordinates(verbose=False)

io = PDBIO()
io.set_structure(structure)
io.save("modified.pdb")
