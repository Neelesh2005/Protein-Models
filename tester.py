import numpy as np

# Define the coordinates from the PDB sequence for ILE A 13
coordinates = {
    'N': np.array([-3.363, 3.528, -18.128]),
    'CA': np.array([-3.507, 2.361, -17.319]),
    'C': np.array([-2.413, 2.231, -16.213]),
    'O': np.array([-2.812, 1.917, -15.092]),
    'CB': np.array([-3.856, 1.134, -18.197]),
    'CG1': np.array([-4.893, 1.637, -19.223]),
    'CG2': np.array([-4.328, -0.049, -17.31]),
    'CD1': np.array([-5.68, 0.461, -19.898]),
    # Add hydrogen atoms if needed
}

# Rotation angle of +90 degrees for psi
theta_deg = 90
theta_rad = np.radians(theta_deg)

# Step 1: Calculate the rotation axis (CA -> C) unit vector
CA = coordinates['CA']
C = coordinates['C']
axis = C - CA
axis = axis / np.linalg.norm(axis)

# Step 2: Create the rotation matrix for +90 degrees
u_x, u_y, u_z = axis
cos_theta = np.cos(theta_rad)
sin_theta = np.sin(theta_rad)
R = np.array([
    [cos_theta + u_x**2 * (1 - cos_theta), u_x * u_y * (1 - cos_theta) - u_z * sin_theta, u_x * u_z * (1 - cos_theta) + u_y * sin_theta],
    [u_y * u_x * (1 - cos_theta) + u_z * sin_theta, cos_theta + u_y**2 * (1 - cos_theta), u_y * u_z * (1 - cos_theta) - u_x * sin_theta],
    [u_z * u_x * (1 - cos_theta) - u_y * sin_theta, u_z * u_y * (1 - cos_theta) + u_x * sin_theta, cos_theta + u_z**2 * (1 - cos_theta)]
])

# Step 3: Apply the rotation matrix to downstream atoms
rotated_coordinates = {}
for atom, coord in coordinates.items():
    if atom in ['CA', 'C','N']:
        # Keep CA and C unchanged
        rotated_coordinates[atom] = coord
    elif coord is not coordinates['N']:  # Apply to downstream atoms (O, CB, CG1, etc.)
        # Translate to the origin at CA
        translated_coord = coord - CA
        # Rotate the translated coordinates
        rotated_translated_coord = R @ translated_coord
        # Translate back to original reference frame
        rotated_coordinates[atom] = rotated_translated_coord + CA

# Output the new coordinates after rotation
print("New coordinates after +90° rotation around the CA-C axis:")
for atom, coord in rotated_coordinates.items():
    print(f"{atom}: {coord}")