import numpy as np
import random
import matplotlib.pyplot as plt

class RubiksCube:
    def __init__(self):
        self.state = np.zeros((3, 3, 3))  # 0 = folded
        self.solved_state = np.zeros((3, 3, 3))  # Reference solved state
    
    def scramble(self):
        self.state = np.random.randint(0, 2, (3, 3, 3))  # Random 0s (folded) and 1s (unfolded)

    def is_solved(self):
        return np.array_equal(self.state, self.solved_state)

    def make_move(self):
        x, y, z = np.random.randint(0, 3, 3)
        self.state[x, y, z] = not self.state[x, y, z]
    
    def energy(self):
        return np.sum(self.state)  # Energy is proportional to the number of unfolded (1s)
    
    def folding_probability(self, temp):
        energy_diff = self.energy()
        return np.exp(-energy_diff / (0.001987 * temp))

def monte_carlo_simulation(cube, num_steps, temperature):
    energy_over_time = []
    for step in range(num_steps):
        cube.make_move()
        fold_prob = cube.folding_probability(temperature)
        
        if random.random() > fold_prob:
            cube.make_move()  # Undo the move if not accepted
        
        energy_over_time.append(cube.energy())
    
    return energy_over_time

# Initialize a Rubik's Cube and scramble it
cube = RubiksCube()
cube.scramble()

# Set simulation parameters
num_steps = 10000
temperature = 600

# Run the Monte Carlo simulation
energy_over_time = monte_carlo_simulation(cube, num_steps, temperature)
if cube.is_solved():
    print("The protein has folded successfully!")
else:
    print("The protein is still partially unfolded.")

# Plot the energy over time to visualize the folding process
plt.plot(energy_over_time)
plt.title('Monte Carlo Simulation of Protein Folding (Rubik\'s Cube Model)')
plt.xlabel('Monte Carlo Steps')
plt.ylabel('Energy (Number of Unfolded States)')
plt.show()
