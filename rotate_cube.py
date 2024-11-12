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
        with open("rotated_coordinates_LEU_4.txt", "a") as file:
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
ATOM    154  N   LEU A  21      -3.609   1.050  -6.006  1.00  0.00           N  
ATOM    155  CA  LEU A  21      -2.631   0.506  -5.161  1.00  0.00           C  
ATOM    156  C   LEU A  21      -2.174   1.424  -4.005  1.00  0.00           C  
ATOM    157  O   LEU A  21      -2.114   0.950  -2.883  1.00  0.00           O  
ATOM    158  CB  LEU A  21      -1.457  -0.290  -5.846  1.00  0.00           C  
ATOM    159  CG  LEU A  21      -1.750  -1.677  -6.376  1.00  0.00           C  
ATOM    160  CD1 LEU A  21      -2.229  -2.585  -5.225  1.00  0.00           C  
ATOM    161  CD2 LEU A  21      -2.741  -1.738  -7.518  1.00  0.00           C  
ATOM    162  H   LEU A  21      -3.287   1.105  -6.948  1.00  0.00           H  
ATOM    163  HA  LEU A  21      -3.189  -0.263  -4.648  1.00  0.00           H  
ATOM    164  HB2 LEU A  21      -0.999   0.346  -6.634  1.00  0.00           H  
ATOM    165  HB3 LEU A  21      -0.632  -0.399  -5.109  1.00  0.00           H  
ATOM    166  HG  LEU A  21      -0.773  -2.097  -6.699  1.00  0.00           H  
ATOM    167 HD11 LEU A  21      -3.309  -2.335  -5.145  1.00  0.00           H  
ATOM    168 HD12 LEU A  21      -1.721  -2.335  -4.269  1.00  0.00           H  
ATOM    169 HD13 LEU A  21      -2.102  -3.672  -5.420  1.00  0.00           H  
ATOM    170 HD21 LEU A  21      -2.829  -2.810  -7.794  1.00  0.00           H  
ATOM    171 HD22 LEU A  21      -2.377  -1.339  -8.490  1.00  0.00           H  
ATOM    172 HD23 LEU A  21      -3.649  -1.177  -7.209  1.00  0.00           H     """  
    
    seq_data = pdb_to_coords.pdb_to_coordinates(pdb_data=pdb_data)
    
    rotator = RotateCube(seq_data=seq_data, amino_cube=cube)
    
    print("Initial cube state:")
    cube.display_cube()
    
    print("\nInitial atom coordinates:")
    for atom, coord in seq_data['coordinates'].items():
        print(f"{atom}: {coord}")
    
    print("\nPerforming rotations:")
    rotator.constant_rotation(7000)

    print("Number of times the sequence was rotated:", rotator.rotation_count)