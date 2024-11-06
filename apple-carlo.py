import yfinance as yf
import numpy as np
import matplotlib.pyplot as plt

stock_symbol = 'CLVT'  # Apple Inc.
data = yf.download(stock_symbol, start='2010-01-01', end='2024-10-28').dropna()

data['Returns'] = data['Adj Close'].pct_change()

mean_return = data['Returns'].mean()
volatility = data['Returns'].std()

num_simulations = 1000  # Number of simulated price paths
num_days = 252  # Simulate for 252 trading days (1 year)

last_price = data['Adj Close'].iloc[-1]

simulated_prices = np.zeros(num_simulations)

for i in range(num_simulations):
    prices = [last_price]  # Start with the last known price
    for _ in range(num_days):
        daily_return = np.random.normal(mean_return, volatility)
        next_price = prices[-1] * (1 + daily_return)
        prices.append(next_price)
    simulated_prices[i] = prices[-1]  # Save the final price

    if i < 10:  # Plot only the first 10 simulations
        plt.plot(prices)

# Plotting settings
plt.title(f"Monte Carlo Simulation for {stock_symbol} Stock Price ({num_simulations} Simulations)")
plt.xlabel("Days")
plt.ylabel("Price")
plt.show()

mean_simulated_price = np.mean(simulated_prices)
max_simulated_price = np.max(simulated_prices)
min_simulated_price = np.min(simulated_prices)

print(f"Mean simulated price after {num_days} days: ${mean_simulated_price:.2f}")
print(f"Max simulated price: ${max_simulated_price:.2f}")
print(f"Min simulated price: ${min_simulated_price:.2f}")
