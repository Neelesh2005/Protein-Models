import numpy as np
import matplotlib.pyplot as plt

protein_length = 50  
num_steps = 1000 
temperature = 300  
protein_state = np.random.choice([0, 1], size=protein_length)

def folding_probability(state, temp):
    energy_diff = -1 if state == 0 else 1  
    return np.exp(-energy_diff / (0.001987 * temp))  
folded_protein_counts = []  
for step in range(num_steps):
    amino_acid = np.random.randint(0, protein_length)

    current_state = protein_state[amino_acid]

    fold_prob = folding_probability(current_state, temperature)

    if np.random.rand() < fold_prob:
        protein_state[amino_acid] = 1 if current_state == 0 else 0

    folded_protein_counts.append(np.sum(protein_state))
folded_count = np.sum(protein_state)
print(f"Final number of folded amino acids: {folded_count} / {protein_length}")
plt.plot(folded_protein_counts)
plt.title("Monte Carlo Simulation of Protein Folding")
plt.xlabel("Monte Carlo Steps")
plt.ylabel("Number of Folded Amino Acids")
plt.show()

folded_count = np.sum(protein_state)
print(f"Final number of folded amino acids: {folded_count} / {protein_length}")
