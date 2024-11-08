import yfinance as yf
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Download data for Indian stocks (Reliance, TCS, and NIFTY 50 ETF)
tickers = ['RELIANCE.NS', 'TCS.NS', 'ICICINIFTY.NS']  # Reliance, TCS, and ICICI NIFTY 50 ETF as baseline
data = yf.download(tickers, start='2023-01-01', end='2023-12-31')['Adj Close']

# Initial portfolio weights
initial_weights = np.array([0.4, 0.3, 0.3])  # Initial weights: 40% Reliance, 30% TCS, 30% NIFTY ETF

# Define momentum-based rotation function
def rotate_allocation_momentum(weights, signal_strength):
    angle = np.radians(45 * signal_strength)  # Rotate based on signal strength, capped at ±45 degrees
    rotation_matrix = np.array([
        [np.cos(angle), -np.sin(angle), 0],
        [np.sin(angle), np.cos(angle), 0],
        [0, 0, 1]
    ])
    rotated_weights = rotation_matrix @ weights
    return np.clip(rotated_weights / np.sum(rotated_weights), 0, 1)  # Normalize

# Calculate a momentum signal (e.g., based on Reliance's 20-day and 50-day moving averages)
short_ma = data['RELIANCE.NS'].rolling(window=20).mean()
long_ma = data['RELIANCE.NS'].rolling(window=50).mean()
momentum_signal = np.where(short_ma > long_ma, 1, -1)  # 1 for positive momentum, -1 for negative

# Apply the dynamic rotation strategy based on momentum signal
rotated_portfolio_values = []
weights = initial_weights.copy()

for i in range(len(data)):
    signal_strength = momentum_signal[i] if not np.isnan(momentum_signal[i]) else 0  # Default to 0 if NaN
    weights = rotate_allocation_momentum(weights, signal_strength)  # Rotate allocation based on momentum
    portfolio_value = (data.iloc[i] * weights).sum()  # Calculate portfolio value
    rotated_portfolio_values.append(portfolio_value)

# Convert portfolio values to a DataFrame for analysis
print(weights)
rotated_portfolio_df = pd.DataFrame({'Portfolio Value': rotated_portfolio_values}, index=data.index)

# Plot the predicted portfolio performance over time
plt.figure(figsize=(10, 6))
plt.plot(rotated_portfolio_df, label='Portfolio Value with Momentum-Based Rotation')
plt.title('Future Portfolio Prediction with Adaptive Weight Rotation (Indian Stocks)')
plt.xlabel('Date')
plt.ylabel('Portfolio Value (INR)')
plt.legend()
plt.show()
