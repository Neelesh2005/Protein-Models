import random
import numpy as np

class ProteinSegment:
    def __init__(self, sequence):
        self.sequence = sequence  
        self.states = self.initialize_states()

    def initialize_states(self):
        # Initialize with random torsion angles (φ, ψ, χ) for each amino acid
        return {i: {"phi": random.uniform(-180, 180), "psi": random.uniform(-180, 180)} 
                for i in range(len(self.sequence))}

    def calculate_energy(self):
        # Placeholder function to calculate the energy of the current conformation
        # Use something like Lennard-Jones potential or a molecular force field
        energy = 0
        for i in range(len(self.sequence) - 1):
            # Dummy energy function based on dihedral angles
            energy += (self.states[i]["phi"] - self.states[i+1]["phi"])**2 + (self.states[i]["psi"] - self.states[i+1]["psi"])**2
        return energy

    def apply_move(self):
        # Randomly change dihedral angles (φ, ψ) for one amino acid (like a Rubik's move)
        i = random.randint(0, len(self.sequence) - 1)
        old_state = self.states[i].copy()
        self.states[i]["phi"] = random.uniform(-180, 180)
        self.states[i]["psi"] = random.uniform(-180, 180)
        return i, old_state

    def revert_move(self, i, old_state):
        # Revert to the previous state if the move is not accepted
        self.states[i] = old_state

    def metropolis_criteria(self, delta_energy, temperature):
        # Metropolis criterion for accepting or rejecting a move
        if delta_energy < 0:
            return True
        else:
            probability = np.exp(-delta_energy / temperature)
            return random.random() < probability

    def run_simulation(self, steps, temperature):
        current_energy = self.calculate_energy()

        for step in range(steps):
            i, old_state = self.apply_move()  # Random move
            new_energy = self.calculate_energy()
            delta_energy = new_energy - current_energy

            if self.metropolis_criteria(delta_energy, temperature):
                current_energy = new_energy  # Accept the new state
            else:
                self.revert_move(i, old_state)  # Revert to old state

            if step % 100 == 0:
                print(f"Step {step}, Energy: {current_energy}")

# Example of running the simulation
protein_sequence = [
    "Ala",  # Alanine
    "Arg",  # Arginine
    "Asn",  # Asparagine
    "Asp",  # Aspartic Acid
    "Cys",  # Cysteine
    "Glu",  # Glutamic Acid
    "Gln",  # Glutamine
    "Gly",  # Glycine
    "His",  # Histidine
    "Ile",  # Isoleucine
    "Leu",  # Leucine
    "Lys",  # Lysine
    "Met",  # Methionine
    "Phe",  # Phenylalanine
    "Pro",  # Proline
    "Ser",  # Serine
    "Thr",  # Threonine
    "Trp",  # Tryptophan
    "Tyr",  # Tyrosine
    "Val",  # Valine
    "Ala",  # Alanine (repeated)
    "Gly",  # Glycine (repeated)
    "Leu",  # Leucine (repeated)
    "Ile",  # Isoleucine (repeated)
    "Ser",  # Serine (repeated)
    "Pro",  # Proline (repeated)
    "Asp"   # Aspartic Acid (repeated)
]

protein = ProteinSegment(protein_sequence)
print(protein.states)
# protein.run_simulation(steps=1000, temperature=300)
