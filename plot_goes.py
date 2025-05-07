import pandas as pd
import matplotlib.pyplot as plt

# Load CSV file
csv_filename = "goes_xray_lightcurve_20240601.csv"
df = pd.read_csv(csv_filename, parse_dates=['Time'], index_col='Time')

# Plot the two GOES channels (short: 0.5–4 Å, long: 1–8 Å)
plt.figure(figsize=(10, 6))
plt.plot(df.index, df['xrsa'], label='0.5–4 Å (short)', color='red')
plt.plot(df.index, df['xrsb'], label='1–8 Å (long)', color='blue')

plt.yscale('log')  # GOES flux is typically plotted on a log scale
plt.xlabel('Time (UTC)')
plt.ylabel('Flux (W/m²)')
plt.title('GOES Soft X-ray Light Curve')
plt.legend()
plt.grid(True)
plt.tight_layout()

plt.show()


