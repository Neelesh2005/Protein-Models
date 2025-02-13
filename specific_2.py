from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.internal_coords import IC_Chain
from io import StringIO

# Properly formatted PDB data with closed triple quotes
pdb_data = '''ATOM      1  N   LEU A  12      -1.987   6.167 -19.323  1.00  0.00           N  
ATOM      2  CA  LEU A  12      -3.221   6.090 -18.481  1.00  0.00           C  
ATOM      3  C   LEU A  12      -3.427   4.823 -17.632  1.00  0.00           C  
ATOM      4  O   LEU A  12      -3.703   4.912 -16.384  1.00  0.00           O  
ATOM      5  CB  LEU A  12      -4.470   6.447 -19.303  1.00  0.00           C  
ATOM      6  CG  LEU A  12      -5.875   6.650 -18.567  1.00  0.00           C  
ATOM      7  CD1 LEU A  12      -5.756   7.830 -17.553  1.00  0.00           C  
ATOM      8  CD2 LEU A  12      -7.010   6.753 -19.546  1.00  0.00           C  
ATOM      9  H   LEU A  12      -1.705   7.138 -19.564  1.00  0.00           H  
ATOM     10  HA  LEU A  12      -3.062   6.746 -17.637  1.00  0.00           H  
ATOM     11  HB2 LEU A  12      -4.259   7.372 -19.882  1.00  0.00           H  
ATOM     12  HB3 LEU A  12      -4.635   5.642 -20.050  1.00  0.00           H  
ATOM     13  HG  LEU A  12      -6.137   5.723 -18.013  1.00  0.00           H  
ATOM     14 HD11 LEU A  12      -5.034   7.636 -16.732  1.00  0.00           H  
ATOM     15 HD12 LEU A  12      -6.727   8.035 -17.053  1.00  0.00           H  
ATOM     16 HD13 LEU A  12      -5.410   8.748 -18.074  1.00  0.00           H  
ATOM     17 HD21 LEU A  12      -7.182   5.785 -20.063  1.00  0.00           H  
ATOM     18 HD22 LEU A  12      -6.690   7.504 -20.300  1.00  0.00           H  
ATOM     19 HD23 LEU A  12      -7.975   7.068 -19.095  1.00  0.00           H  
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
'''

def rotate_phi(residue_id: int, chain_id: str = "A"):
    # Read structure from string
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("input", r"E:\protein-model\8U1T_correct.pdb")
    
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                ic_chain = IC_Chain(chain)
                ic_chain.atom_to_internal_coordinates()
                
                # Find target residue
                target_res = chain[(' ', residue_id, ' ')]
                if target_res.internal_coord:
                    current_phi = target_res.internal_coord.get_angle("phi")
                    if current_phi is not None:
                        new_phi = (current_phi + 90) % 360
                        target_res.internal_coord.set_angle("phi", new_phi)
                
                # Update coordinates
                ic_chain.internal_to_atom_coordinates()
    
    # Save modified structure
    io = PDBIO()
    io.set_structure(structure)
    io.save(f"modified_phi_{residue_id}.pdb")

# Rotate phi angle for ILE A13 (residue ID 13)
rotate_phi(13)
