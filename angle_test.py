from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.vectors import calc_dihedral

# Load the protein structure
parser = PDBParser(QUIET=True)
structure = parser.get_structure("8U1T", r"E:\Downloads\8u1t (1).pdb")

for model in structure:
    for chain in model:
        print(f"Chain {chain.id}, length {len(chain)}")
        residues = list(chain)
        print(f"From {residues[0].get_resname()} to {residues[-1].get_resname()}")
        
        # Calculate phi and psi angles for each residue
        for i, residue in enumerate(residues[1:-1], 1):
            try:
                phi = calc_dihedral(
                    residues[i-1]["C"].get_vector(),
                    residue["N"].get_vector(),
                    residue["CA"].get_vector(),
                    residue["C"].get_vector()
                )
                psi = calc_dihedral(
                    residue["N"].get_vector(),
                    residue["CA"].get_vector(),
                    residue["C"].get_vector(),
                    residues[i+1]["N"].get_vector()
                )
                print(f"{residue.get_resname()}: phi = {phi * 57.2958:.2f}, psi = {psi * 57.2958:.2f}")
            except KeyError:
                # Skip if atoms are missing
                print(f"{residue.get_resname()}: phi/psi angles not defined")
