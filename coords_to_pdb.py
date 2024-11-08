import numpy as np
import pdb_to_coords
pdb_data = """ATOM     20  N   ILE A  13      -3.363   3.528 -18.128  1.00  0.00           N  
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
def coordinates_to_pdb(coordinates, pdb_info):
    pdb_lines = []
    
    for info in pdb_info:
        atom_number = info["atom_number"]
        atom_name = info["atom_name"]
        residue_name = info["residue_name"]
        chain_id = info["chain_id"]
        residue_seq = info["residue_seq"]
        element = info["element"]
        
        # Retrieve the coordinates for the current atom
        coord = coordinates[atom_name]
        x, y, z = coord[0], coord[1], coord[2]
        
        # Format to PDB line
        pdb_line = f"ATOM  {atom_number:5d} {atom_name:<4} {residue_name} {chain_id} {residue_seq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {element}"
        pdb_lines.append(pdb_line)
    
    return "\n".join(pdb_lines)

# Example usage
