import numpy as np
import random
import single_amino_2
import rot_90_phi
import rot_90_psi
import pdb_to_coords
import coords_to_pdb
import json

class RotateCube:
    def __init__(self, seq_data: dict, amino_cube: single_amino_2.AminoCube):
        self.amino_cube = amino_cube
        self.cube = amino_cube.get_cube()
        self.sequence_cords = amino_cube.cords
        self.sequence_name = amino_cube.sequence_name
        self.atom_Cords = seq_data.get('coordinates')
        self.atom_Info = seq_data.get('info')
        self.rotation_count = 0
        self.structure_data = {}

    def select_random_dimension(self):
        return random.choice(["x", "y", "z"])

    def select_random_layer(self):
        return random.randint(0, 2)
    
    def direction(self):
        return random.choice([1, -1])

    def save_coordinates_to_json(self):
        pdb_sequence = coords_to_pdb.coordinates_to_pdb(self.atom_Cords, self.atom_Info)

    # Create a new structure number based on the current rotation count
        structure_number = str(self.rotation_count)

    # Load existing data if the file exists, else create a new dictionary
        try:
            with open("rotated_coordinates.json", "r") as file:
                data = json.load(file)
        except FileNotFoundError:
            data = {}

    # Store the PDB sequence as a multi-line string
        data[structure_number] = pdb_sequence.strip()

    # Save the updated data back to the JSON file with formatting
        with open("test.json", "w") as file:
            json.dump(data, file, indent=4)
    def constant_rotation(self, number):
        for _ in range(number):
            print("Step:", _ + 1)
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
            print('Rotated and saving to JSON file')
            self.rotation_count += 1
            self.save_coordinates_to_json()

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
ATOM    352  N   ILE A  33      -2.707   1.602  12.056  1.00  0.00           N  
ATOM    353  CA  ILE A  33      -2.742   2.579  13.188  1.00  0.00           C  
ATOM    354  C   ILE A  33      -4.063   2.635  13.928  1.00  0.00           C  
ATOM    355  O   ILE A  33      -4.025   2.955  15.117  1.00  0.00           O  
ATOM    356  CB  ILE A  33      -2.293   3.991  12.842  1.00  0.00           C  
ATOM    357  CG1 ILE A  33      -3.074   4.662  11.674  1.00  0.00           C  
ATOM    358  CG2 ILE A  33      -0.798   3.951  12.389  1.00  0.00           C  
ATOM    359  CD1 ILE A  33      -4.367   5.359  11.962  1.00  0.00           C  
ATOM    360  H   ILE A  33      -2.727   1.898  11.104  1.00  0.00           H  
ATOM    361  HA  ILE A  33      -2.060   2.286  13.973  1.00  0.00           H  
ATOM    362  HB  ILE A  33      -2.202   4.651  13.731  1.00  0.00           H  
ATOM    363 HG12 ILE A  33      -2.324   5.342  11.216  1.00  0.00           H  
ATOM    364 HG13 ILE A  33      -3.277   3.892  10.900  1.00  0.00           H  
ATOM    365 HG21 ILE A  33      -0.151   3.481  13.159  1.00  0.00           H  
ATOM    366 HG22 ILE A  33      -0.367   4.962  12.224  1.00  0.00           H  
ATOM    367 HG23 ILE A  33      -0.743   3.256  11.524  1.00  0.00           H  
ATOM    368 HD11 ILE A  33      -4.955   5.419  11.022  1.00  0.00           H  
ATOM    369 HD12 ILE A  33      -5.026   4.805  12.665  1.00  0.00           H  
ATOM    370 HD13 ILE A  33      -4.156   6.380  12.344  1.00  0.00           H      """
    seq_data = pdb_to_coords.pdb_to_coordinates(pdb_data=pdb_data)
    
    rotator = RotateCube(seq_data=seq_data, amino_cube=cube)
    
    print("Initial cube state:")
    cube.display_cube()
    
    print("\nInitial atom coordinates:")
    for atom, coord in seq_data['coordinates'].items():
        print(f"{atom}: {coord}")
    
    print("\nPerforming rotations:")
    rotator.constant_rotation(1000)

    print("Number of times the sequence was rotated:", rotator.rotation_count)
