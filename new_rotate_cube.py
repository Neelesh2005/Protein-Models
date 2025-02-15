import numpy as np
import random
import new_cube
import new_cube_helix_preserved_rot_90_phi as helix_preserved_rot_90_phi
import new_cube_helix_preserved_rot_90_psi as helix_preserved_rot_90_psi
from Bio.PDB import PDBParser, PDBIO
class RotateCube:
    def __init__(self,pdb_file,chain_id, amino_cube: new_cube.AminoCube):
        self.amino_cube = amino_cube
        self.cube = amino_cube.get_cube()
        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure("helix", pdb_file)
        self.chain_id = chain_id
        self.sequence_name = amino_cube.sequence_name
        self.rotation_count = 0
        self.res_num = amino_cube.res_num
        self.unique_structures = set()
        
    def select_random_dimension(self):
        return random.choice(["x", "y", "z"])

    def select_random_layer(self):
        return random.randint(0, 2)
    
    def direction(self):
        return random.choice([1, -1])

    def save_coordinates_to_file(self):
        io = PDBIO()
        io.set_structure(self.structure)
        output_file = f"helix_rotated_psi_{self.sequence_name}.pdb"
        io.save(output_file)
        print(f"Successfully rotated structure saved to {output_file}")
        
    def constant_rotation(self,number):
        for _ in range(number):
            print("Step:", _+1)
            print("rotating now")
            self.perform_random_rotation()
        # print(f"Number of unique structures: {len(self.unique_structures)}")
        self.save_coordinates_to_file()
    # def no_file_constant_rotation(self, number):
    #     for _ in range(number):
            
    #         self.no_file_rotation()
    #         atom_positions = tuple(tuple(coord) for coord in self.atom_Cords.values())
    #         self.unique_structures.add(atom_positions)
        
    def rotate_layer(self, axis, layer, direction):
        rotated = False 

        if axis == "x":
            self.cube[layer, :, :] = np.rot90(self.cube[layer, :, :], k=direction)
            
        elif axis == "y":
            self.cube[:, layer, :] = np.rot90(self.cube[:, layer, :], k=direction)
            if self.sequence_name in self.cube[:, layer, :]:
                self.rotation_count +=1
                prev_struct = self.structure
                print("rotated by phi angle")
                # print("Before rotation:", [(atom.coord) for atom in prev_struct.get_atoms()])

                self.structure = helix_preserved_rot_90_phi.rotate_phi_preserve_helix(prev_struct,self.chain_id,self.res_num,direction*np.random.randint(90,180))
                # print("After rotation:", [(atom.coord) for atom in self.structure.get_atoms()])
                rotated = True
                
        elif axis == "z":
            self.cube[:, :, layer] = np.rot90(self.cube[:, :, layer], k=direction)
            if self.sequence_name in self.cube[:, :, layer]:
                prev_struct = self.structure
                self.rotation_count +=1
                # print("Before rotation:", [(atom.coord) for atom in prev_struct.get_atoms()])
                print("rotated by psi angle")
                self.structure = helix_preserved_rot_90_psi.rotate_psi_preserve_helix(prev_struct,self.chain_id,self.res_num,direction*np.random.randint(90,180))
                # print("After rotation:", [(atom.coord) for atom in self.structure.get_atoms()])
                rotated = True

        # if rotated:
        #     print('Rotated')
        #     print('saving to file')
        #     self.rotation_count += 1
        #     self.save_coordinates_to_file()

        # print(f"Rotated layer {layer} along {axis}-axis with direction {direction}")
    # def no_file_rotation(self):
    #     axis = self.select_random_dimension()
    #     layer = self.select_random_layer()
    #     direction = self.direction()
    #     def no_file_rotate_layer(axis, layer, direction):
    #         rotated = False

    #         if axis == "x":
    #             self.cube[layer, :, :] = np.rot90(self.cube[layer, :, :], k=direction)
                
    #         elif axis == "y":
    #             self.cube[:, layer, :] = np.rot90(self.cube[:, layer, :], k=direction)
    #             if self.sequence_name in self.cube[:, layer, :]:
    #                 self.atom_Cords = rot_90_psi.rotate_coordinates(self.atom_Cords, direction * 90)
    #                 rotated = True
                    
    #         elif axis == "z":
    #             self.cube[:, :, layer] = np.rot90(self.cube[:, :, layer], k=direction)
    #             if self.sequence_name in self.cube[:, :, layer]:
    #                 self.atom_Cords = rot_90_phi.rotate_coordinates(self.atom_Cords, direction * 90)
    #                 rotated = True

    #         if rotated:
    #             self.rotation_count += 1
                
    #     no_file_rotate_layer(axis, layer, direction)
    def perform_random_rotation(self):
        axis = self.select_random_dimension()
        layer = self.select_random_layer()
        direction = self.direction()
        print(axis,layer,direction)
        self.rotate_layer(axis, layer, direction)

if __name__ == "__main__":
    cube = new_cube.AminoCube(r"e:\\protein-model\\8U1T_correct.pdb","A")
    cube.place_at_random()
    cube.display_cube()
    rotator = RotateCube(r"e:\\protein-model\\8U1T_correct.pdb","A", cube)
    # rotator.perform_random_rotation()
    rotator.constant_rotation(100)
    print(rotator.rotation_count)
    # seq_data = pdb_to_coords.pdb_to_coordinates(pdb_data=pdb_data)
    
    # rotator = RotateCube(seq_data=seq_data, amino_cube=cube)
    
    # print("Initial cube state:")
    # cube.display_cube()
    
    # print("\nInitial atom coordinates:")
    # for atom, coord in seq_data['coordinates'].items():
    #     print(f"{atom}: {coord}")
    
    # print("\nPerforming rotations:")
    # rotator.constant_rotation(100)

    # print("Number of times the sequence was rotated:", rotator.rotation_count)