import numpy as np
from itertools import combinations, groupby

class BondLineStructure:
    def __init__(self, pdb_data):
        self.atoms = []
        self.bonds = []
        self.parse_pdb_data(pdb_data)

    def parse_pdb_data(self, pdb_data):
        """Parse PDB data and store atom information."""
        for line in pdb_data.strip().split('\n'):
            if line.startswith('ATOM'):
                atom = {
                    'id': int(line[6:11].strip()),
                    'name': line[12:16].strip(),
                    'element': line[76:78].strip(),
                    'coords': np.array([
                        float(line[30:38].strip()),
                        float(line[38:46].strip()),
                        float(line[46:54].strip())
                    ])
                }
                self.atoms.append(atom)

    def calculate_distance(self, atom1, atom2):
        """Calculate distance between two atoms."""
        return np.linalg.norm(atom1['coords'] - atom2['coords'])

    def identify_bonds(self, max_bond_length=1.6):
        """Identify bonds between atoms based on distance."""
        for atom1, atom2 in combinations(self.atoms, 2):
            if self.calculate_distance(atom1, atom2) <= max_bond_length:
                self.bonds.append((atom1, atom2))

    def generate_bond_line_structure(self):
        """Generate a simplified bond-line structure."""
        self.identify_bonds()

        # Group atoms by residue (useful for clearer output)
        residue_groups = groupby(self.atoms, key=lambda x: x['name'])

        structure = []
        visited = set()  # To avoid duplicate bonds

        for atom1, atom2 in self.bonds:
            if (atom1['id'], atom2['id']) not in visited and (atom2['id'], atom1['id']) not in visited:
                # Create a bond-line notation
                bond = f"{atom1['element']}-{atom2['element']}"
                structure.append(bond)
                visited.add((atom1['id'], atom2['id']))

        return " + ".join(structure)

# Example usage
pdb_data = """
ATOM      1  N   LEU A  12      -1.987   6.167 -19.323  1.00  0.00           N  
ATOM      2  CA  LEU A  12      -3.221   6.090 -18.481  1.00  0.00           C  
ATOM      3  C   LEU A  12      -3.427   4.823 -17.632  1.00  0.00           C  
ATOM      4  O   LEU A  12      -3.703   4.912 -16.384  1.00  0.00           O  
ATOM      5  CB  LEU A  12      -4.470   6.447 -19.303  1.00  0.00           C  
ATOM      6  CG  LEU A  12      -5.875   6.650 -18.567  1.00  0.00           C  
ATOM      7  CD1 LEU A  12      -5.756   7.830 -17.553  1.00  0.00           C  
ATOM      8  CD2 LEU A  12      -7.010   6.753 -19.546  1.00  0.00           C  
ATOM      9  H   LEU A  12      -1.705   7.138 -19.564  1.00  0.00           H  
ATOM     10  HA  LEU A  12      -3.062   6.746 -17.637  1.00  0.00           H  
"""

analyzer = BondLineStructure(pdb_data)
result = analyzer.generate_bond_line_structure()
print("Bond-Line Structure:", result)
