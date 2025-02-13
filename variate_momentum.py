# import yfinance as yf
# import pandas as pd
# import numpy as np

# # Fetch historical stock data
# def get_stock_data(ticker, start_date, end_date):
#     stock_data = yf.download(ticker, start=start_date, end=end_date)
#     return stock_data

# # Calculate momentum
# def calculate_momentum(data, period=14):
#     data['Momentum'] = data['Close'] - data['Close'].shift(period)
#     return data

# # Calculate binomial variate
# def calculate_binomial_variate(data, period=14):
#     data['Up'] = np.where(data['Momentum'] > 0, 1, 0)
#     data['Down'] = np.where(data['Momentum'] <= 0, 1, 0)
#     data['Binomial'] = (data['Up'] - data['Down']).rolling(window=period).sum()
#     return data

# # Generate trading signals based on binomial variate
# def generate_signals(data, threshold=5):
#     data['Signal'] = 0
#     data.loc[data['Binomial'] > threshold, 'Signal'] = 1  # Buy signal
#     data.loc[data['Binomial'] < -threshold, 'Signal'] = -1  # Sell signal
#     return data

# # Backtest the strategy
# def backtest_strategy(data):
#     data['Returns'] = data['Close'].pct_change()
#     data['Strategy_Returns'] = data['Returns'] * data['Signal'].shift(1)
#     data['Cumulative_Returns'] = (1 + data['Strategy_Returns']).cumprod()
#     return data

# # Parameters
# ticker = 'AAPL'
# start_date = '2020-01-01'
# end_date = '2025-01-01'
# momentum_period = 14
# binomial_threshold = 5

# # Run the strategy
# stock_data = get_stock_data(ticker, start_date, end_date)
# stock_data = calculate_momentum(stock_data, momentum_period)
# stock_data = calculate_binomial_variate(stock_data, momentum_period)
# stock_data = generate_signals(stock_data, binomial_threshold)
# stock_data = backtest_strategy(stock_data)

# # Display the results
# print(stock_data[['Close', 'Momentum', 'Binomial', 'Signal', 'Cumulative_Returns']].tail())

# # Plot the results
# import matplotlib.pyplot as plt

# plt.figure(figsize=(12, 6))
# plt.plot(stock_data['Cumulative_Returns'], label='Strategy Returns')
# plt.plot((1 + stock_data['Returns']).cumprod(), label='Buy and Hold Returns')
# plt.legend()
# plt.show()

import pandas as pd
import yfinance as yf
from arch import arch_model
import numpy as np
# Fetch historical stock data
ticker = 'AAPL'
start_date = '2020-01-01'
end_date = '2025-01-01'
stock_data = yf.download(ticker, start=start_date, end=end_date)

# Calculate log returns
stock_data['Log_Returns'] = np.log(stock_data['Close'] / stock_data['Close'].shift(1))
returns = stock_data['Log_Returns'].dropna()

# Fit GARCH(1, 1) model
model = arch_model(returns, vol='Garch', p=1, q=1)
model_fit = model.fit(disp='off')
print(model_fit.summary())

# Forecast volatility
volatility_forecast = model_fit.forecast(horizon=5)
print(volatility_forecast.variance[-1:])