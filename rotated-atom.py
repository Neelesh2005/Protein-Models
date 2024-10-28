import numpy as np

def rotation_matrix(axis, theta):
    axis = axis / np.linalg.norm(axis)
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    return np.array([
        [a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d + a*c)],
        [2*(b*c + a*d), a*a + c*c - b*b - d*d, 2*(c*d - a*b)],
        [2*(b*d - a*c), 2*(c*d + a*b), a*a + d*d - b*b - c*c]
    ])

# Function to rotate an atom around a bond
def rotate_atom(atom, bond_start, bond_end, angle):
    # Translate so the bond_start is at the origin
    translated_atom = atom - bond_start
    translated_bond_end = bond_end - bond_start

    # Rotation matrix
    rotation_mat = rotation_matrix(translated_bond_end, angle)

    # Rotate the atom
    rotated_atom = np.dot(rotation_mat, translated_atom)

    # Translate back
    return rotated_atom + bond_start

# Atom coordinates based on input data
atoms = {
    "O": np.array([-3.703, 4.912, -16.384]),
    "CB": np.array([-4.470, 6.447, -19.303]),
    "CG": np.array([-5.875, 6.650, -18.567]),
    "CD1": np.array([-5.756, 7.830, -17.553]),
    "CD2": np.array([-7.010, 6.753, -19.546]),
    "H": np.array([-1.705, 7.138, -19.564]),
    "HA": np.array([-3.062, 6.746, -17.637]),
    "HB2": np.array([-4.259, 7.372, -19.882]),
    "HB3": np.array([-4.635, 5.642, -20.050]),
    "HG": np.array([-6.137, 5.723, -18.013]),
    "HD11": np.array([-5.034, 7.636, -16.732])
}

# Define the bonds between atoms (key atoms for side chain rotation)
bonds = [
    ("CB", "CG"),  # Rotate around CB-CG bond
    ("CG", "CD1"), # Rotate around CG-CD1 bond
    ("CG", "CD2")  # Rotate around CG-CD2 bond
]

# Angle increments (in radians)
angle_increment = np.pi / 6  # 30 degrees

# Function to apply rotations around specified bonds
def generate_conformations(atoms, bonds, angle_increment):
    conformations = []
    for bond in bonds:
        bond_start = atoms[bond[0]]
        bond_end = atoms[bond[1]]
        
        # Rotate by multiples of 30 degrees (12 steps for 360 degrees)
        for i in range(12):
            angle = i * angle_increment
            rotated_atoms = {}
            for atom_name, atom_coord in atoms.items():
                # Only rotate atoms beyond the bond (atoms after bond_end)
                if atom_name != bond[0] and atom_name != bond[1]:
                    rotated_atoms[atom_name] = rotate_atom(atom_coord, bond_start, bond_end, angle)
                else:
                    rotated_atoms[atom_name] = atom_coord  # Keep these atoms fixed

            conformations.append(rotated_atoms)
    return conformations

# Generate all conformations
conformations = generate_conformations(atoms, bonds, angle_increment)

# Print out the generated conformations
for idx, conformation in enumerate(conformations):
    print(f"Conformation {idx + 1}:")
    for atom, coord in conformation.items():
        print(f"{atom}: {coord}")
    print("\n")
