import random
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
class ProteinSegment:
    def __init__(self, sequence):
        self.sequence = sequence
        self.states = self.initialize_states()
        self.positions = self.initialize_positions()

    def initialize_states(self):
        return {i: {"phi": random.uniform(-180, 180), "psi": random.uniform(-180, 180)}
                for i in range(len(self.sequence))}

    def initialize_positions(self):
        positions = {}
        for i in range(len(self.sequence)):
            positions[i] = np.array([i, 0, 0])
        return positions

    def calculate_energy(self):
        energy = 0
        #not accurate use CHARMM fu]orce field
        for i in range(len(self.sequence) - 1):
            distance = np.linalg.norm(self.positions[i] - self.positions[i + 1])
            energy += (1 / distance) ** 12 - 2 * (1 / distance) ** 6
        return energy

    def apply_move(self):
        i = random.randint(0, len(self.sequence) - 1)
        old_state = self.states[i].copy()
        old_position = self.positions[i].copy()

        self.states[i]["phi"] = random.uniform(-180, 180)
        self.states[i]["psi"] = random.uniform(-180, 180)

        new_position = self.calculate_new_position(i)

        if not self.is_valid_move(new_position):
            return i, old_state, old_position, False

        self.positions[i] = new_position
        return i, old_state, old_position, True

    def calculate_new_position(self, i):
        direction = np.array([
            np.cos(np.radians(self.states[i]["phi"])),
            np.sin(np.radians(self.states[i]["phi"])),
            np.sin(np.radians(self.states[i]["psi"]))
        ])

        if i == 0:  # Handle the first amino acid differently
            new_position = self.positions[i] + direction
        else:
            new_position = self.positions[i - 1] + direction  # Use the previous position for others

        return new_position

    def is_valid_move(self, new_position):
        for j in range(len(self.positions)):
            if np.array_equal(self.positions[j], new_position):
                return False
        return True

    def revert_move(self, i, old_state, old_position):
        self.states[i] = old_state
        self.positions[i] = old_position

    def energy_criteria(self, delta_energy, temperature):
        if delta_energy < 0:
            return True
        else:
            probability = np.exp(-delta_energy / temperature)
            return random.random() < probability

    def run_simulation(self, steps, temperature):
        current_energy = self.calculate_energy()
        energies = [current_energy]  # Store the initial energy

        for step in range(steps):
            i, old_state, old_position, valid_move = self.apply_move()

            if not valid_move:
                continue

            new_energy = self.calculate_energy()
            delta_energy = new_energy - current_energy

            if self.energy_criteria(delta_energy, temperature):
                current_energy = new_energy
            else:
                self.revert_move(i, old_state, old_position)

            energies.append(current_energy)  # Store energy at every step

            if step % 10 == 0:
                print(f"Step {step}, Energy: {current_energy},phi {self.states[i]["phi"]}, psi: {self.states[i]["psi"]}, Number: {i}")

        plt.figure(figsize=(10, 6))
        plt.plot(energies, label='Energy', color='b')
        plt.xlabel('Steps')
        plt.ylabel('Energy')
        plt.ylim(-50, 200)

        plt.title('Energy vs. Simulation Steps')
        plt.legend()
        plt.grid(True)
        plt.show()

        return energies  # Return the array of energies if needed

# Example of running the simulation
protein_sequence = [
    "Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", "Ile", 
    "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", 
    "Ala", "Gly", "Leu", "Ile", "Ser", "Pro", "Asp"
]

protein = ProteinSegment(protein_sequence)

energies = protein.run_simulation(steps=10000, temperature=300)
np_energy = np.array(energies)
test = pd.DataFrame(np_energy)
print(float(np.min(np_energy[100:9900])) ,test[100:9900].idxmin())