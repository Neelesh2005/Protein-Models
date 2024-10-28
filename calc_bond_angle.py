import numpy as np
from itertools import combinations

class ProteinStructureAnalyzer:
    def __init__(self, pdb_data):
        self.atoms = []
        self.parse_pdb_data(pdb_data)
        
    def parse_pdb_data(self, pdb_data):
        for line in pdb_data.strip().split('\n'):
            if line.startswith('ATOM'):
                atom = {
                    'id': int(line[6:11].strip()),
                    'name': line[12:16].strip(),
                    'residue': line[17:20].strip(),
                    'chain': line[21].strip(),
                    'residue_num': int(line[22:26].strip()),
                    'coords': np.array([
                        float(line[30:38].strip()),
                        float(line[38:46].strip()),
                        float(line[46:54].strip())
                    ]),
                    'element': line[76:78].strip()
                }
                self.atoms.append(atom)

    def calculate_distance(self, atom1, atom2):
        return np.linalg.norm(atom1['coords'] - atom2['coords'])

    def calculate_bond_angle(self, atom1, atom2, atom3):
        vector1 = atom1['coords'] - atom2['coords']
        vector2 = atom3['coords'] - atom2['coords']
        
        cos_angle = np.dot(vector1, vector2) / (np.linalg.norm(vector1) * np.linalg.norm(vector2))
        cos_angle = np.clip(cos_angle, -1.0, 1.0)  # Handle floating-point precision issues
        angle = np.arccos(cos_angle)
        
        return np.degrees(angle)

    def identify_bonded_atoms(self, max_bond_length=100000.0):
        bonds = []
        for atom1, atom2 in combinations(self.atoms, 2):
            distance = self.calculate_distance(atom1, atom2)
            if distance <= max_bond_length:
                bonds.append((atom1, atom2))
        return bonds

    def analyze_bond_angles(self):
        bonds = self.identify_bonded_atoms()
        
        # Create a dictionary to store bonded atoms for each atom
        bonded_atoms = {atom['id']: [] for atom in self.atoms}
        for atom1, atom2 in bonds:
            bonded_atoms[atom1['id']].append(atom2)
            bonded_atoms[atom2['id']].append(atom1)

        # Calculate angles for all bonded triplets
        angles = []
        for center_atom in self.atoms:
            bonded_to_center = bonded_atoms[center_atom['id']]
            for atom1, atom2 in combinations(bonded_to_center, 2):
                angle = self.calculate_bond_angle(atom1, center_atom, atom2)
                angles.append({
                    'atom1': atom1,
                    'center': center_atom,
                    'atom3': atom2,
                    'angle': angle
                })
        
        return angles

    def print_analysis(self):
        print("Protein Structure Analysis")
        
        # Print atom information
        print("Atom Information:")
        for atom in self.atoms:
            print(f"Atom {atom['id']}: {atom['name']} ({atom['element']}) in {atom['residue']} {atom['residue_num']}")
            print(f"  Coordinates: ({atom['coords'][0]:.3f}, {atom['coords'][1]:.3f}, {atom['coords'][2]:.3f})")

        # Print bonds
        print("\nBond Information:")
        bonds = self.identify_bonded_atoms()
        if bonds:
            for atom1, atom2 in bonds:
                distance = self.calculate_distance(atom1, atom2)
                print(f"Bond: {atom1['name']}{atom1['id']} - {atom2['name']}{atom2['id']}")
                print(f"  Distance: {distance:.3f} Å")
        else:
            print("No bonds found.")

        # Print angles
        print("\nBond Angles:")
        angles = self.analyze_bond_angles()
        if angles:
            for angle_info in angles:
                atom1 = angle_info['atom1']
                center = angle_info['center']
                atom3 = angle_info['atom3']
                print(f"Angle: {atom1['name']}{atom1['id']} - {center['name']}{center['id']} - {atom3['name']}{atom3['id']}")
                print(f"  {angle_info['angle']:.2f} degrees")
        else:
            print("No angles found.")

# Example usage
if __name__ == "__main__":
    pdb_data = """
    ATOM      1  N   LEU A  12      -1.987   6.167 -19.323  1.00  0.00           N  
    ATOM      2  CA  LEU A  12      -3.221   6.090 -18.481  1.00  0.00           C  
    ATOM      3  C   LEU A  12      -3.427   4.823 -17.632  1.00  0.00           C  
    ATOM      4  O   LEU A  12      -3.703   4.912 -16.384  1.00  0.00           O  
    ATOM      5  CB  LEU A  12      -4.470   6.447 -19.303  1.00  0.00           C  
    ATOM      6  CG  LEU A  12      -5.875   6.650 -18.567  1.00  0.00           C  
    ATOM      7  CD1 LEU A  12      -5.756   7.830 -17.553  1.00  0.00           C  
    ATOM      8  CD2 LEU A  12      -7.010   6.753 -19.546  1.00  0.00           C  
    """
    
    analyzer = ProteinStructureAnalyzer(pdb_data)
    analyzer.print_analysis()
