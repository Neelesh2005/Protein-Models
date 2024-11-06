import numpy as np
from typing import Dict, List, Tuple

class AtomCoordinate:
    def __init__(self, x: float, y: float, z: float):
        self.coord = np.array([x, y, z])
    
    def get_vector(self):
        return self.coord

class SingleResidueCalculator:
    def __init__(self, pdb_text: str):
        """Initialize with PDB text format for a single residue."""
        self.atoms = self._parse_pdb_text(pdb_text)
        
    def _parse_pdb_text(self, pdb_text: str) -> Dict[str, AtomCoordinate]:
        """Parse PDB text to extract atom coordinates."""
        atoms = {}
        for line in pdb_text.split('\n'):
            if line.startswith('ATOM'):
                atom_name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms[atom_name] = AtomCoordinate(x, y, z)
        return atoms
    
    def _calculate_dihedral(self, p1: np.ndarray, p2: np.ndarray, 
                          p3: np.ndarray, p4: np.ndarray) -> float:
        """Calculate dihedral angle between four points."""
        b1 = p2 - p1
        b2 = p3 - p2
        b3 = p4 - p3
        
        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)
        
        m1 = np.cross(n1, b2/np.linalg.norm(b2))
        
        x = np.dot(n1, n2)
        y = np.dot(m1, n2)
        
        angle = np.arctan2(y, x)
        return np.degrees(angle)

    def calculate_chi_angles(self) -> Dict[str, float]:
        """Calculate chi angles for the residue."""
        # Chi angle definitions for LEU
        chi_angles = {}
        
        try:
            # Chi1: N-CA-CB-CG
            if all(atom in self.atoms for atom in ['N', 'CA', 'CB', 'CG']):
                chi1 = self._calculate_dihedral(
                    self.atoms['N'].get_vector(),
                    self.atoms['CA'].get_vector(),
                    self.atoms['CB'].get_vector(),
                    self.atoms['CG'].get_vector()
                )
                chi_angles['chi1'] = chi1
            
            # Chi2: CA-CB-CG-CD1/CD2
            if all(atom in self.atoms for atom in ['CA', 'CB', 'CG', 'CD1']):
                chi2 = self._calculate_dihedral(
                    self.atoms['CA'].get_vector(),
                    self.atoms['CB'].get_vector(),
                    self.atoms['CG'].get_vector(),
                    self.atoms['CD1'].get_vector()
                )
                chi_angles['chi2'] = chi2
            
        except KeyError as e:
            print(f"Warning: Missing atom {e} for angle calculation")
        
        return chi_angles

def format_angles(angles: Dict[str, float], precision: int = 2) -> str:
    """Format angles dictionary into a readable string."""
    output = ["Calculated angles for LEU:"]
    for angle_type, value in angles.items():
        if value is not None:
            output.append(f"{angle_type}: {value:.{precision}f}°")
        else:
            output.append(f"{angle_type}: N/A")
    return "\n".join(output)

# Example usage
if __name__ == "__main__":
    # Your PDB text
    pdb_text = """MODEL        1                                                                  
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
ATOM     11  HB2 LEU A  12      -4.259   7.372 -19.882  1.00  0.00           H  
ATOM     12  HB3 LEU A  12      -4.635   5.642 -20.050  1.00  0.00           H  
ATOM     13  HG  LEU A  12      -6.137   5.723 -18.013  1.00  0.00           H  
ATOM     14 HD11 LEU A  12      -5.034   7.636 -16.732  1.00  0.00           H  
ATOM     15 HD12 LEU A  12      -6.727   8.035 -17.053  1.00  0.00           H  
ATOM     16 HD13 LEU A  12      -5.410   8.748 -18.074  1.00  0.00           H  
ATOM     17 HD21 LEU A  12      -7.182   5.785 -20.063  1.00  0.00           H  
ATOM     18 HD22 LEU A  12      -6.690   7.504 -20.300  1.00  0.00           H  
ATOM     19 HD23 LEU A  12      -7.975   7.068 -19.095  1.00  0.00           H"""
    
    calculator = SingleResidueCalculator(pdb_text)
    angles = calculator.calculate_chi_angles()
    print(format_angles(angles))