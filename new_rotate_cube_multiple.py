import numpy as np
import random
import new_cube_multiple as new_cube
import new_cube_helix_preserved_rot_90_phi as helix_preserved_rot_90_phi
import new_cube_helix_preserved_rot_90_psi as helix_preserved_rot_90_psi
from Bio.PDB import PDBParser, PDBIO

class RotateCube:
    def __init__(self, pdb_file, chain_id, amino_cube: new_cube.AminoCube):
        self.amino_cube = amino_cube
        self.cube = amino_cube.get_cube()
        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure("helix", pdb_file)
        self.chain_id = chain_id
        self.sequence_names = amino_cube.sequence_names
        self.rotation_count = 0
        self.res_num = amino_cube.res_nums
        self.unique_structures = set()
        self.residue_map = dict(zip(self.sequence_names, self.res_num))
        
    def select_random_dimension(self):
        return random.choice(["x", "y", "z"])

    def select_random_layer(self):
        return random.randint(0, 2)
    
    def direction(self):
        return random.choice([1, -1])

    def save_coordinates_to_file(self):
        io = PDBIO()
        io.set_structure(self.structure)
        output_file = f"helix_rotated_{random.randint(100,999)}.pdb"
        io.save(output_file)
        print(f"Successfully rotated structure saved to {output_file}")
        
    def constant_rotation(self, number):
        for _ in range(number):
            print("Step:", _+1)
            print("Rotating now")
            self.perform_random_rotation()
        self.save_coordinates_to_file()
    
    def rotate_layer(self, axis, layer, direction):
        rotated_residues = set()

        if axis == "x":
            self.cube[layer, :, :] = np.rot90(self.cube[layer, :, :], k=direction)

        elif axis == "y":
            self.cube[:, layer, :] = np.rot90(self.cube[:, layer, :], k=direction)
            layer_residues = set(self.cube[:, layer, :].flatten())

            for res in layer_residues:
                if res in self.residue_map and res not in rotated_residues:
                    self.rotation_count += 1
                    prev_struct = self.structure
                    print(f"Rotated {res} by phi angle")
                    self.structure = helix_preserved_rot_90_phi.rotate_phi_preserve_helix(
                        prev_struct, self.chain_id, self.residue_map[res], direction * np.random.randint(90, 180)
                    )
                    rotated_residues.add(res)

        elif axis == "z":
            self.cube[:, :, layer] = np.rot90(self.cube[:, :, layer], k=direction)
            layer_residues = set(self.cube[:, :, layer].flatten())

            for res in layer_residues:
                if res in self.residue_map and res not in rotated_residues:
                    self.rotation_count += 1
                    prev_struct = self.structure
                    print(f"Rotated {res} by psi angle")
                    self.structure = helix_preserved_rot_90_psi.rotate_psi_preserve_helix(
                        prev_struct, self.chain_id, self.residue_map[res], direction * np.random.randint(90, 180)
                    )
                    rotated_residues.add(res)

    def perform_random_rotation(self):
        axis = self.select_random_dimension()
        layer = self.select_random_layer()
        direction = self.direction()
        print(f"Rotating axis={axis}, layer={layer}, direction={direction}")
        self.rotate_layer(axis, layer, direction)

if __name__ == "__main__":
    cube = new_cube.AminoCube(r"e:\\protein-model\\8U1T_correct.pdb", "A")
    cube.place_at_random()
    cube.display_cube()
    rotator = RotateCube(r"e:\\protein-model\\8U1T_correct.pdb", "A", cube)
    rotator.constant_rotation(100)
    print("Total rotations:", rotator.rotation_count)
