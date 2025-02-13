import numpy as np



def pdb_to_coordinates(pdb_data):
    coordinates = {}
    pdb_info = []

    for line in pdb_data.strip().splitlines():
        parts = line.split()
        
        atom_number = int(parts[1])
        atom_name = parts[2]
        residue_name = parts[3]
        chain_id = parts[4]
        residue_seq = int(parts[5])
        x, y, z = map(float, parts[6:9])
        element = parts[-1]  # The last part is usually the element symbol

        # Store coordinates
        coordinates[atom_name] = np.array([x, y, z])

        # Collect necessary PDB information to reconstruct the sequence
        pdb_info.append({
            "atom_number": atom_number,
            "atom_name": atom_name,
            "residue_name": residue_name,
            "chain_id": chain_id,
            "residue_seq": residue_seq,
            "element": element
        })

    return {"coordinates":coordinates, "info":pdb_info}

# Convert PDB data to coordinates and pdb info
