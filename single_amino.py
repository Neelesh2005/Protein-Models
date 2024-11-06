import numpy as np
import random

class AminoCube:
    def __init__(self, number, sequence="ALA"):
        # Initialize an empty 3x3x3 matrix and set attributes
        self.number = number
        self.sequence = sequence
        self.reset_cube()

    def reset_cube(self):
        # Reset the cube to an empty state
        self.cube = np.random.randint(1, 101, size=(3, 3, 3))

    def place_at_edge_corner(self):
        # Clear the cube before placing
        self.reset_cube()
        # Randomly place at one of the 8 corners of the cube
        x = random.choice([0, 2])  # Front or back depth layer
        y = random.choice([0, 2])  # Bottom or top height layer
        z = random.choice([0, 2])  # Left or right width layer
        self.cords = [x, y, z]
        self.cube[x, y, z] = self.sequence
        return self.cords

    def place_at_edge_center(self):
        # Clear the cube before placing
        self.reset_cube()
        # Place at one of the 12 edge centers (1 in one dimension, 0 or 2 in the others)
        x, y, z = random.choice([(1, 0, 0), (1, 0, 2), (1, 2, 0), (1, 2, 2),
                                 (0, 1, 0), (0, 1, 2), (2, 1, 0), (2, 1, 2),
                                 (0, 0, 1), (2, 0, 1), (0, 2, 1), (2, 2, 1)])
        self.cords = [x, y, z]
        self.cube[x, y, z] = self.sequence
        return self.cords

    def place_at_center(self):
        # Clear the cube before placing
        self.reset_cube()
        # Place within the central structure (center layer in depth, avoiding edges and corners)
        x = 1  # Middle depth layer
        y, z = random.choice([(0, 1), (2, 1), (1, 0), (1, 2)])  # Avoiding corners and edges
        self.cords = [x, y, z]
        self.cube[x, y, z] = self.sequence
        return self.cords

    def place_at_random(self):
        # Clear the cube before placing
        self.reset_cube()
        # Place anywhere in the cube randomly
        x, y, z = random.randint(0, 2), random.randint(0, 2), random.randint(0, 2)
        self.cords = [x, y, z]
        self.cube[x, y, z] = self.sequence
        return self.cords

    def place_at_coords(self, x, y, z):
        # Clear the cube before placing
        self.reset_cube()
        # Place at specified coordinates with validation
        if 0 <= x <= 2 and 0 <= y <= 2 and 0 <= z <= 2:
            self.cords = [x, y, z]
            self.cube[x, y, z] = self.sequence
            return self.cords
        else:
            raise ValueError("Coordinates must be between 0 and 2 for a 3x3x3 cube.")

    def display_cube(self):
        # Display the cube
        print("3x3x3 Matrix with current placements:")
        print(self.cube)


# Example usage
cube = AminoCube(1, sequence="A")

# Place amino acid at different positions and print the results
cube.place_at_edge_corner()
print("After placing at an edge corner:")
cube.display_cube()

cube.place_at_edge_center()
print("After placing at an edge center:")
cube.display_cube()

cube.place_at_center()
print("After placing at the center structure (central layer, excluding edges):")
cube.display_cube()
