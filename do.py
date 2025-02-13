import numpy as np
import matplotlib.pyplot as plt

# Define price range
S_T = np.linspace(80, 120, 100)  # Future spot price at expiry

# Given futures contract price (K)
K = 100

# Payoff for short futures position: Profit = K - S_T
payoff = K - S_T

# Plot the payoff diagram
plt.figure(figsize=(8, 5))
plt.plot(S_T, payoff, label="Short Futures Payoff", color="red", linewidth=2)
plt.axhline(0, color='black', linewidth=1, linestyle="--")  # Zero profit line
plt.axvline(K, color='gray', linestyle="--", label="Futures Price (K)")

# Labels and title
plt.xlabel("Spot Price at Expiry (S_T)")
plt.ylabel("Profit / Loss")
plt.title("Payoff Diagram for Short Futures Position")
plt.legend()
plt.grid(True)

# Show the graph
plt.show()