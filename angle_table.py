import numpy as np
import pdb_to_coords

def calculate_phi_psi(sequence_dict):
    """
    Calculates phi and psi angles for each sequence in the dictionary.
    
    Args:
        sequence_dict (dict): Dictionary with sequence names as keys and PDB data as values.
    """
    results = {}

    for seq_name, pdb_data in sequence_dict.items():
        print(f"\nCalculating phi and psi angles for sequence: {seq_name}\n")
        
        # Convert PDB data to atomic coordinates
        seq_data = pdb_to_coords.pdb_to_coordinates(pdb_data=pdb_data)

        phi_angles = []
        psi_angles = []

        seq_coordinates = seq_data['coordinates']
        
        keys = list(seq_coordinates.keys())
        
        for i in range(1, len(keys) - 1):  # Avoid first and last residues
            # Extract atomic coordinates
            try:
                prev_c_coords = seq_coordinates[keys[i - 1]]
                n_coords = seq_coordinates[keys[i]]
                ca_coords = seq_coordinates[keys[i]]
                ca_coords = seq_coordinates[keys[i]]
                c_coords = seq_coordinates[keys[i]]
                next_n_coords = seq_coordinates[keys[i + 1]]

                # Calculate phi and psi
                phi = calculate_dihedral(prev_c_coords, n_coords, ca_coords, c_coords)
                psi = calculate_dihedral(n_coords, ca_coords, c_coords, next_n_coords)

                phi_angles.append(phi)
                psi_angles.append(psi)

            except KeyError as e:
                print(f"Missing atom {e} in {seq_name}. Skipping residue.")
            except Exception as e:
                print(f"Error processing residue {keys[i]}: {e}")

        results[seq_name] = {'phi': phi_angles, 'psi': psi_angles}

        # Display statistics
        if phi_angles and psi_angles:
            print(f"Mean phi: {np.mean(phi_angles):.2f}, Variance: {np.var(phi_angles):.2f}")
            print(f"Mean psi: {np.mean(psi_angles):.2f}, Variance: {np.var(psi_angles):.2f}")
        else:
            print(f"No phi/psi angles calculated for {seq_name}.")

    return results


def calculate_dihedral(p1, p2, p3, p4):
    """
    Calculate the dihedral angle between four points using vector math.
    
    Args:
        p1, p2, p3, p4: 3D coordinate tuples (x, y, z)

    Returns:
        Dihedral angle in degrees.
    """
    p1, p2, p3, p4 = np.array(p1), np.array(p2), np.array(p3), np.array(p4)

    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3

    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    n1 /= np.linalg.norm(n1)
    n2 /= np.linalg.norm(n2)

    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)

    return np.degrees(np.arctan2(y, x))


if __name__ == "__main__":
    sequence_dict = {
        "Sequence_1": """
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
ATOM    370 HD13 ILE A  33      -4.156   6.380  12.344  1.00  0.00           H 
        """}
    phi_psi_results = calculate_phi_psi(sequence_dict)

    # Print results
    for seq_name, angles in phi_psi_results.items():
        print(f"\n{seq_name} Phi Angles: {angles['phi']}")
        print(f"{seq_name} Psi Angles: {angles['psi']}")