import numpy as np
import random
import single_amino_2

class RotateCube:
    def __init__(self, amino_cube: single_amino_2.AminoCube):
        self.amino_cube = amino_cube
        self.cube = amino_cube.get_cube()
        self.sequence_cords = amino_cube.cords
        self.sequence_name = amino_cube.sequence_name
        self.phi = amino_cube.phi
        self.psi = amino_cube.psi
        self.chi1 = amino_cube.chi1

    def select_random_dimension(self):
        return random.choice(["x", "y", "z"])

    def select_random_layer(self):
        return random.randint(0, 2)
    
    def direction(self):
        return random.choice([1, -1])

    def rotate_layer(self, axis, layer, direction):
        if axis == "x":
            self.cube[layer, :, :] = np.rot90(self.cube[layer, :, :], k=direction)
            if(self.sequence_name in self.cube[layer, :, :]):
                self.chi1 = (self.chi1 + 90 * direction)
        elif axis == "y":
            self.cube[:, layer, :] = np.rot90(self.cube[:, layer, :], k=direction)
            if(self.sequence_name in self.cube[:, layer, :]):
                self.psi = (self.psi + 90 * direction)
        elif axis == "z":
            self.cube[:, :, layer] = np.rot90(self.cube[:, :, layer], k=direction)
            if(self.sequence_name in self.cube[:, :, layer]):
                self.phi = (self.phi + 90 * direction)

        print(f"Rotated layer {layer} along {axis}-axis with direction {direction}:\n{self.cube}\n, phi = {self.phi}, psi = {self.psi}, chi1 = {self.chi1}")

    def perform_random_rotation(self):
        axis = self.select_random_dimension()
        layer = self.select_random_layer()
        direction = self.direction()
        self.rotate_layer(axis, layer, direction)

if __name__ == "__main__":
    cube = single_amino_2.AminoCube(1, sequence={'name': 'A', 'phi': 180, 'psi': 180, 'chi1': 180})
    cube.place_at_random()  
    
    rotator = RotateCube(cube)
    
    print("Initial cube state:")
    cube.display_cube()
    
    print("\nPerforming rotations:")
    for i in range(3):
        print(f"\nRotation {i+1}:")
        rotator.perform_random_rotation()