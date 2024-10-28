
import random
import math
import matplotlib.pyplot as plt

# Function to estimate π using Monte Carlo simulation
def monte_carlo_pi(num_points):
    inside_circle = 0

    for _ in range(num_points):
        x = random.uniform(0, 1)  # Random x coordinate
        y = random.uniform(0, 1)  # Random y coordinate

        # Check if the point is inside the quarter circle
        if x**2 + y**2 <= 1:
            inside_circle += 1

    # π is approximately 4 times the ratio of points inside the circle to total points
    pi_estimate = 4 * (inside_circle / num_points)
    return pi_estimate

# Simulate with a large number of points
num_points = 10000
pi_estimate = monte_carlo_pi(num_points)
print(f"Estimated value of π after {num_points} points: {pi_estimate}")

# Plot the points for visualization
inside_x = []
inside_y = []
outside_x = []
outside_y = []

for _ in range(num_points):
    x = random.uniform(0, 1)
    y = random.uniform(0, 1)
    if x**2 + y**2 <= 1:
        inside_x.append(x)
        inside_y.append(y)
    else:
        outside_x.append(x)
        outside_y.append(y)

plt.figure(figsize=(6, 6))
plt.scatter(inside_x, inside_y, color='blue', label="Inside Circle", s=1)
plt.scatter(outside_x, outside_y, color='red', label="Outside Circle", s=1)
plt.legend()
plt.title(f"Monte Carlo Simulation with {num_points} Points")
plt.show()
