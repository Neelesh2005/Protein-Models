import numpy as np
import matplotlib.pyplot as plt

# Parameters for the simulation
protein_length = 50  # Number of amino acids in the protein
num_steps = 1000  # Number of Monte Carlo steps
temperature = 300  # Temperature in Kelvin (affects folding probability)

# Initialize the protein in a random unfolded state
protein_state = np.random.choice([0, 1], size=protein_length)

# Folding/unfolding probabilities (simplified Boltzmann factor)
def folding_probability(state, temp):
    energy_diff = -1 if state == 0 else 1  # Energy difference for folded/unfolded states
    return np.exp(-energy_diff / (0.001987 * temp))  # Boltzmann factor (R constant in kcal/mol·K)

# Monte Carlo simulation
folded_protein_counts = []  # Store number of folded amino acids at each step

for step in range(num_steps):
    # Pick a random amino acid to try to fold/unfold
    amino_acid = np.random.randint(0, protein_length)

    # Current state of the amino acid (0 = unfolded, 1 = folded)
    current_state = protein_state[amino_acid]

    # Calculate the folding/unfolding probability
    fold_prob = folding_probability(current_state, temperature)

    # Randomly decide whether to fold/unfold the amino acid
    if np.random.rand() < fold_prob:
        protein_state[amino_acid] = 1 if current_state == 0 else 0

    # Count how many amino acids are folded at this step
    folded_protein_counts.append(np.sum(protein_state))
folded_count = np.sum(protein_state)
print(f"Final number of folded amino acids: {folded_count} / {protein_length}")
# Plot the results of the simulation
plt.plot(folded_protein_counts)
plt.title("Monte Carlo Simulation of Protein Folding")
plt.xlabel("Monte Carlo Steps")
plt.ylabel("Number of Folded Amino Acids")
plt.show()

# Final state analysis
folded_count = np.sum(protein_state)
print(f"Final number of folded amino acids: {folded_count} / {protein_length}")
