import numpy as np

def rotate_coordinates(coordinates:dict, angle_deg):
    angle_rad = np.radians(angle_deg)
    CA = coordinates.get('CA')
    # print(CA)
    C = coordinates.get('C')
    # print(sC)
    axis = C - CA
    axis = axis / np.linalg.norm(axis)

    u_x, u_y, u_z = axis
    cos_theta = np.cos(angle_rad)
    sin_theta = np.sin(angle_rad)
    R = np.array([
        [cos_theta + u_x**2 * (1 - cos_theta), u_x * u_y * (1 - cos_theta) - u_z * sin_theta, u_x * u_z * (1 - cos_theta) + u_y * sin_theta],
        [u_y * u_x * (1 - cos_theta) + u_z * sin_theta, cos_theta + u_y**2 * (1 - cos_theta), u_y * u_z * (1 - cos_theta) - u_x * sin_theta],
        [u_z * u_x * (1 - cos_theta) - u_y * sin_theta, u_z * u_y * (1 - cos_theta) + u_x * sin_theta, cos_theta + u_z**2 * (1 - cos_theta)]
    ])

    rotated_coordinates = dict()
    for atom, coord in coordinates.items():
        if atom in ['CA', 'C', 'N']:
            rotated_coordinates[atom] = coord
        else:
            translated_coord = coord - CA
            rotated_translated_coord = R @ translated_coord
            rotated_coordinates[atom] = rotated_translated_coord + CA

    return rotated_coordinates


# coordinates = {
#     'N': np.array([-3.363, 3.528, -18.128]),
#     'CA': np.array([-3.507, 2.361, -17.319]),
#     'C': np.array([-2.413, 2.231, -16.213]),
#     'O': np.array([-2.812, 1.917, -15.092]),
#     'CB': np.array([-3.856, 1.134, -18.197]),
#     'CG1': np.array([-4.893, 1.637, -19.223]),
#     'CG2': np.array([-4.328, -0.049, -17.31]),
#     'CD1': np.array([-5.68, 0.461, -19.898]),
# }


# rotated_coordinates = rotate_coordinates(coordinates, 90)
# for atom, coord in rotated_coordinates.items():
#     print(f"{atom}: {coord}")
