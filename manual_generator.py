import random

def generate_files(file_dict):
    filenames = list(file_dict.keys())
    num_files = random.randint(1, len(filenames))
    selected_files = random.sample(filenames, num_files)
    generated_files = {}
    for filename in selected_files:
        sequence_count = random.randint(1, file_dict[filename])
        generated_files[filename] = sequence_count

    return generated_files

if __name__ == "__main__":
    file_dict = {
        "rotated_coordinates_ALA_2.txt": 100,
        "rotated_coordinates_ALA.txt": 120,
        "rotated_coordinates_ASN.txt": 150,
        "rotated_coordinates_ILE_2.txt": 80,
        "rotated_coordinates_ILE.txt": 200,
        "rotated_coordinates_LEU_3.txt": 180,
        "rotated_coordinates_LEU_4.txt": 160,
        "rotated_coordinates_LEU_6.txt": 140,
        "rotated_coordinates_LEU.txt": 220,
        "rotated_coordinates_PHE_2.txt": 130,
        "rotated_coordinates_PHE.txt": 170,
        "rotated_coordinates_SER.txt": 110,
        "rotated_coordinates_THR.txt": 190,
        "rotated_coordinates_VAL_2.txt": 210,
        "rotated_coordinates_VAL.txt": 250
    }
    for _ in range(20):
        result = generate_files(file_dict)
        print("Generated Files with Sequence Counts:")
        print(result)
