import numpy as np
import random
from typing import List

class AminoCube:
    def __init__(self, number: int, sequence: dict = {'name':'A','phi':180,'psi':180,'chi1':180}):
        self.number = number
        self.sequence_name = sequence.get('name')
        self.cords = None
        self.phi = sequence.get('phi')
        self.psi = sequence.get('psi')
        self.chi1 = sequence.get('chi1')
        self.reset_cube()

    def reset_cube(self):
        self.cube = np.full((3, 3, 3), None)

    def place_at_edge_corner(self) -> List[int]:
        self.reset_cube()
        # Pick front/back for depth, bottom/top for height, left/right for width
        x = random.choice([0, 2])
        y = random.choice([0, 2])
        z = random.choice([0, 2])
        
        self.cords = [x, y, z]
        self.cube[x, y, z] = self.sequence_name
        return self.cords

    def place_at_edge_center(self) -> List[int]:
        self.reset_cube()
        
        edges = [
            # Front-back edges (middle depth)
            (1, 0, 0), (1, 0, 2), (1, 2, 0), (1, 2, 2),
            # Top-bottom edges (middle height)
            (0, 1, 0), (0, 1, 2), (2, 1, 0), (2, 1, 2),
            # Left-right edges (middle width)
            (0, 0, 1), (2, 0, 1), (0, 2, 1), (2, 2, 1)
        ]
        
        x, y, z = random.choice(edges)
        self.cords = [x, y, z]
        self.cube[x, y, z] = self.sequence_name
        return self.cords

    def place_at_center(self) -> List[int]:
        self.reset_cube()
        
        # Always middle depth, varying height and width
        x = 1
        centers = [(0, 1), (2, 1), (1, 0), (1, 2)]
        y, z = random.choice(centers)
        
        self.cords = [x, y, z]
        self.cube[x, y, z] = self.sequence_name
        return self.cords

    def place_at_random(self) -> List[int]:
        self.reset_cube()
        x = random.randint(0, 2)
        y = random.randint(0, 2)
        z = random.randint(0, 2)
        
        self.cords = [x, y, z]
        self.cube[x, y, z] = self.sequence_name
        return self.cords

    def place_at_coords(self, x: int, y: int, z: int) -> List[int]:
        if not (0 <= x <= 2 and 0 <= y <= 2 and 0 <= z <= 2):
            raise ValueError("Coordinates must be between 0 and 2")
            
        self.reset_cube()
        self.cords = [x, y, z]
        self.cube[x, y, z] = self.sequence
        return self.cords

    def display_cube(self):
        if self.cords:
            print(f"\nPosition: d={self.cords[0]}, h={self.cords[1]}, w={self.cords[2]}")
        
        print("\nCube state:")
        print(self.cube)
        return self.cube  # Added return statement

    def get_cube(self):  # Added new method
        return self.cube

# Usage example
if __name__ == "__main__":
    cube = AminoCube(1, {'name':'A','phi':180,'psi':180,'chi1':180})
    
    print("Corner placement:")
    cube.place_at_edge_corner()
    cube.display_cube()
    
    print("\nEdge center placement:")
    cube.place_at_edge_center()
    cube.display_cube()
    
    print("\nCenter placement:")
    cube.place_at_center()
    cube.display_cube()