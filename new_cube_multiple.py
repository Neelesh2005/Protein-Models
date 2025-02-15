import numpy as np
import random
from typing import List
from Bio.PDB import PDBParser
class AminoCube:
    def __init__(self,pdb_file, chain_id):
        
        self.chain_id = chain_id
        self.structure = PDBParser().get_structure('helix', pdb_file)    
        self.model = self.structure[0]
        chain = self.model[chain_id]
        self.res_list = list(chain.get_residues())
        self.reset_cube()
        self.length = len(self.res_list)
        self.indexes = np.random.randint(0,self.length,size=3)
        self.res_nums = [self.res_list[index].get_id()[1] for index in self.indexes]
        self.sequence_names = [self.res_list[index].get_resname() for index in self.indexes]
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
    def rand_coords(self):
        return random.randint(0, 2), random.randint(0, 2), random.randint(0, 2)
    def place_at_random(self) -> List[int]:
        self.reset_cube()
        prev_cords = None
        for i in range(len(self.sequence_names)):
            x, y, z = self.rand_coords()
            if (x, y, z) == prev_cords:
                continue
            prev_cords = (x, y, z)
            self.cords = [x, y, z]
            self.cube[x, y, z] = self.sequence_names[i]
        

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
   
    cube = AminoCube(r"e:\\protein-model\\8U1T_correct.pdb","A")
    cube.place_at_random()
    cube.display_cube()
    print(cube.sequence_names)
    print(cube.res_nums)
    # print("Corner placement:")
    # cube.place_at_edge_corner()
    # cube.display_cube()
    
    # print("\nEdge center placement:")
    # cube.place_at_edge_center()
    # cube.display_cube()
    
    # print("\nCenter placement:")
    # cube.place_at_center()
    # cube.display_cube()