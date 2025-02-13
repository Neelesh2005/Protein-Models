import single_amino_2
import pdb_to_coords
from rotate_cube import RotateCube
import matplotlib.pyplot as plt
import numpy as np

def run_simulations(iterations, rotation_steps):
    pdb_data = """
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
    """
    seq_data = pdb_to_coords.pdb_to_coordinates(pdb_data=pdb_data)
    
    unique_counts = []

    for i in range(iterations):
        cube = single_amino_2.AminoCube(1, sequence={'name': 'A', 'phi': 180, 'psi': 180, 'chi1': 180})
        cube.place_at_random()
        
        rotator = RotateCube(seq_data=seq_data, amino_cube=cube)
        rotator.no_file_constant_rotation(rotation_steps)

        unique_count = len(rotator.unique_structures)
        unique_counts.append(unique_count)
        
        print(f"Iteration {i + 1}: {unique_count} unique structures")

    # Calculate statistics
    mean_count = np.mean(unique_counts)
    variance_count = np.var(unique_counts)

    print("\nSummary of Results:")
    for idx, count in enumerate(unique_counts):
        print(f"Iteration {idx + 1}: {count} unique structures")
    print(f"\nMean of unique structures: {mean_count}")
    print(f"Variance of unique structures: {variance_count}")

    # Plot the results
    plt.figure(figsize=(10, 6))
    plt.plot(range(1, iterations + 1), unique_counts, marker='o', linestyle='-', color='b', label='Unique Structures')
    plt.axhline(y=mean_count, color='r', linestyle='--', label=f"Mean ({mean_count:.2f})")
    plt.title("Unique Structures vs. Iterations")
    plt.xlabel("Iteration Number")
    plt.ylabel("Number of Unique Structures")
    plt.legend()
    plt.grid()
    plt.show()

if __name__ == "__main__":
    iterations = 1000  # Number of simulations to run
    rotation_steps = 7000  # Number of rotation steps per simulation
    run_simulations(iterations, rotation_steps)
