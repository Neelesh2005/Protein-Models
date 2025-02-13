import numpy as np
import random

# Define the PDB data as a constant since it wasn't being passed correctly
PDB_DATA = """
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
ATOM     19 HD23 LEU A  12      -7.975   7.068 -19.095  1.00  0.00           H  
ATOM     20  N   ILE A  13      -3.363   3.528 -18.128  1.00  0.00           N  
ATOM     21  CA  ILE A  13      -3.507   2.361 -17.319  1.00  0.00           C  
ATOM     22  C   ILE A  13      -2.413   2.231 -16.213  1.00  0.00           C  
ATOM     23  O   ILE A  13      -2.812   1.917 -15.092  1.00  0.00           O  
ATOM     24  CB  ILE A  13      -3.856   1.134 -18.197  1.00  0.00           C  
ATOM     25  CG1 ILE A  13      -4.893   1.637 -19.223  1.00  0.00           C  
ATOM     26  CG2 ILE A  13      -4.328  -0.049 -17.310  1.00  0.00           C  
ATOM     27  CD1 ILE A  13      -5.680   0.461 -19.898  1.00  0.00           C  
ATOM     28  H   ILE A  13      -3.122   3.308 -19.070  1.00  0.00           H  
ATOM     29  HA  ILE A  13      -4.343   2.626 -16.688  1.00  0.00           H  
ATOM     30  HB  ILE A  13      -2.937   0.808 -18.728  1.00  0.00           H  
ATOM     31 HG12 ILE A  13      -5.595   2.363 -18.759  1.00  0.00           H  
ATOM     32 HG13 ILE A  13      -4.460   2.172 -20.095  1.00  0.00           H  
ATOM     33 HG21 ILE A  13      -3.665  -0.294 -16.453  1.00  0.00           H  
ATOM     34 HG22 ILE A  13      -4.681  -0.883 -17.954  1.00  0.00           H  
ATOM     35 HG23 ILE A  13      -5.335   0.139 -16.881  1.00  0.00           H  
ATOM     36 HD11 ILE A  13      -4.925  -0.256 -20.287  1.00  0.00           H  
ATOM     37 HD12 ILE A  13      -6.287  -0.046 -19.117  1.00  0.00           H  
ATOM     38 HD13 ILE A  13      -6.288   0.616 -20.815  1.00  0.00           H  
ATOM     39  N   VAL A  14      -1.131   2.537 -16.456  1.00  0.00           N  
ATOM     40  CA  VAL A  14      -0.117   2.522 -15.317  1.00  0.00           C  
ATOM     41  C   VAL A  14      -0.426   3.491 -14.235  1.00  0.00           C  
ATOM     42  O   VAL A  14      -0.288   3.210 -13.050  1.00  0.00           O  
ATOM     43  CB  VAL A  14       1.357   2.606 -15.744  1.00  0.00           C  
ATOM     44  CG1 VAL A  14       2.481   2.903 -14.752  1.00  0.00           C  
ATOM     45  CG2 VAL A  14       1.683   1.384 -16.613  1.00  0.00           C  
ATOM     46  H   VAL A  14      -0.752   2.693 -17.365  1.00  0.00           H  
ATOM     47  HA  VAL A  14      -0.175   1.575 -14.802  1.00  0.00           H  
ATOM     48  HB  VAL A  14       1.431   3.495 -16.406  1.00  0.00           H  
ATOM     49 HG11 VAL A  14       3.475   2.810 -15.239  1.00  0.00           H  
ATOM     50 HG12 VAL A  14       2.392   2.178 -13.915  1.00  0.00           H  
ATOM     51 HG13 VAL A  14       2.337   3.922 -14.336  1.00  0.00           H  
ATOM     52 HG21 VAL A  14       2.758   1.414 -16.892  1.00  0.00           H  
ATOM     53 HG22 VAL A  14       1.130   1.262 -17.569  1.00  0.00           H  
ATOM     54 HG23 VAL A  14       1.491   0.411 -16.114  1.00  0.00           H  
"""

class PDBMatrix:
    def __init__(self):
        self.matrix = np.empty((6, 3, 3), dtype=dict)
        self.parse_pdb_data(PDB_DATA)
    
    def parse_pdb_data(self, pdb_data):
        # Split the data by lines and filter empty lines
        pdb_lines = [line for line in pdb_data.split('\n') if line.startswith('ATOM')]
        
        # Fill the matrix
        for i in range(6):
            for j in range(3):
                for k in range(3):
                    idx = i * 9 + j * 3 + k
                    if idx < len(pdb_lines):
                        line = pdb_lines[idx]
                        self.matrix[i][j][k] = {
                            "atom": line[12:16].strip(),
                            "residue": line[17:20].strip(),
                            "coords": np.array([
                                float(line[30:38].strip()),
                                float(line[38:46].strip()),
                                float(line[46:54].strip())
                            ]),
                            "element": line[76:78].strip()
                        }
    
    def rotate_coords(self, coords, axis, angle_degrees):
        """Rotate coordinates around specified axis by given angle."""
        angle = np.radians(angle_degrees)
        c = np.cos(angle)
        s = np.sin(angle)
        
        if axis == 'x':
            rotation_matrix = np.array([
                [1, 0, 0],
                [0, c, -s],
                [0, s, c]
            ])
        elif axis == 'y':
            rotation_matrix = np.array([
                [c, 0, s],
                [0, 1, 0],
                [-s, 0, c]
            ])
        else:  # z-axis
            rotation_matrix = np.array([
                [c, -s, 0],
                [s, c, 0],
                [0, 0, 1]
            ])
            
        return np.dot(rotation_matrix, coords)
    
    def rotate_layer_z(self, layer):
        """Rotate coordinates in a z-layer by 90 degrees."""
        if 0 <= layer < 6:
            for i in range(3):
                for j in range(3):
                    if self.matrix[layer][i][j] is not None:
                        self.matrix[layer][i][j]['coords'] = self.rotate_coords(
                            self.matrix[layer][i][j]['coords'], 'z', 90
                        )
    
    def rotate_layer_x(self, layer):
        """Rotate coordinates in an x-layer by 90 degrees."""
        if 0 <= layer < 3:
            for i in range(6):
                for j in range(3):
                    if self.matrix[i][j][layer] is not None:
                        self.matrix[i][j][layer]['coords'] = self.rotate_coords(
                            self.matrix[i][j][layer]['coords'], 'x', 90
                        )
    
    def rotate_layer_y(self, layer):
        """Rotate coordinates in a y-layer by 90 degrees."""
        if 0 <= layer < 3:
            for i in range(6):
                for k in range(3):
                    if self.matrix[i][layer][k] is not None:
                        self.matrix[i][layer][k]['coords'] = self.rotate_coords(
                            self.matrix[i][layer][k]['coords'], 'y', 90
                        )
    
    def random_rotation(self):
        """Perform a random rotation if probabilities are approximately equal."""
        p1 = random.random()
        p2 = random.random()
        
        if abs(p1 - p2) < 0.05:
            axis = random.choice(['x', 'y', 'z'])
            if axis == 'x':
                layer = random.randint(0, 2)
                self.rotate_layer_x(layer)
                return f"Rotated coordinates in layer {layer} along x-axis"
            elif axis == 'y':
                layer = random.randint(0, 2)
                self.rotate_layer_y(layer)
                return f"Rotated coordinates in layer {layer} along y-axis"
            else:
                layer = random.randint(0, 5)
                self.rotate_layer_z(layer)
                return f"Rotated coordinates in layer {layer} along z-axis"
        return "No rotation performed"
    
    def print_matrix(self):
        """Print the current state of the matrix."""
        for i in range(6):
            print(f"\nFace {i+1}:")
            for j in range(3):
                row = []
                for k in range(3):
                    atom = self.matrix[i][j][k]
                    if atom is not None:
                        coords = atom['coords']
                        row.append(f"{atom['atom']}({coords[0]:.2f},{coords[1]:.2f},{coords[2]:.2f})")
                    else:
                        row.append("None")
                print(row)

# Example usage
if __name__ == "__main__":
    # Create the PDB matrix
    pdb_matrix = PDBMatrix()
    
    # Print initial state
    print("Initial matrix state:")
    pdb_matrix.print_matrix()
    
    # Perform random rotation
    result = pdb_matrix.random_rotation()
    print(f"\n{result}")
    
    # Print final state
    print("\nFinal matrix state:")
    pdb_matrix.print_matrix()