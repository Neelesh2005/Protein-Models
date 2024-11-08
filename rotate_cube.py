import numpy as np
import random
import single_amino_2
import rot_90_phi
import rot_90_psi
import pdb_to_coords
import coords_to_pdb
class RotateCube:
    def __init__(self, seq_data: dict, amino_cube: single_amino_2.AminoCube):
        self.amino_cube = amino_cube
        self.cube = amino_cube.get_cube()
        self.sequence_cords = amino_cube.cords
        self.sequence_name = amino_cube.sequence_name
        self.atom_Cords = seq_data.get('coordinates')
        self.atom_Info = seq_data.get('info')
        self.rotation_count = 0

    def select_random_dimension(self):
        return random.choice(["x", "y", "z"])

    def select_random_layer(self):
        return random.randint(0, 2)
    
    def direction(self):
        return random.choice([1, -1])

    def save_coordinates_to_file(self):
        with open("rotated_coordinates.txt", "a") as file:
            pdb_sequence = coords_to_pdb.coordinates_to_pdb(self.atom_Cords, self.atom_Info)
    
            file.write("Updated PDB sequence after rotation:\n")
            file.write(pdb_sequence)
            file.write("\n\n")
    def constant_rotation(self,number):
        for _ in range(number):
            print("Step:", _+1)
            self.perform_random_rotation()
    def rotate_layer(self, axis, layer, direction):
        rotated = False 

        if axis == "x":
            self.cube[layer, :, :] = np.rot90(self.cube[layer, :, :], k=direction)
            
        elif axis == "y":
            self.cube[:, layer, :] = np.rot90(self.cube[:, layer, :], k=direction)
            if self.sequence_name in self.cube[:, layer, :]:
                self.atom_Cords = rot_90_psi.rotate_coordinates(self.atom_Cords, direction * 90)
                rotated = True
                
        elif axis == "z":
            self.cube[:, :, layer] = np.rot90(self.cube[:, :, layer], k=direction)
            if self.sequence_name in self.cube[:, :, layer]:
                self.atom_Cords = rot_90_phi.rotate_coordinates(self.atom_Cords, direction * 90)
                rotated = True

        if rotated:
            print('Rotated')
            print('saving to file')
            self.rotation_count += 1
            self.save_coordinates_to_file()

        print(f"Rotated layer {layer} along {axis}-axis with direction {direction}")
        
    def perform_random_rotation(self):
        axis = self.select_random_dimension()
        layer = self.select_random_layer()
        direction = self.direction()
        
        self.rotate_layer(axis, layer, direction)

if __name__ == "__main__":
    cube = single_amino_2.AminoCube(1, sequence={'name': 'A', 'phi': 180, 'psi': 180, 'chi1': 180})
    cube.place_at_random()  
    pdb_data = """
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
ATOM     38 HD13 ILE A  13      -6.288   0.616 -20.815  1.00  0.00           H  """  
    
    seq_data = pdb_to_coords.pdb_to_coordinates(pdb_data=pdb_data)
    
    rotator = RotateCube(seq_data=seq_data, amino_cube=cube)
    
    print("Initial cube state:")
    cube.display_cube()
    
    print("\nInitial atom coordinates:")
    for atom, coord in seq_data['coordinates'].items():
        print(f"{atom}: {coord}")
    
    print("\nPerforming rotations:")
    rotator.constant_rotation(10000)

    print("Number of times the sequence was rotated:", rotator.rotation_count)