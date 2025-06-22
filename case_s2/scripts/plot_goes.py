import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()

# Load CSV file
csv_filename = "/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case6_Oct09/goes_xray_lightcurve_20241009.csv"
df = pd.read_csv(csv_filename, parse_dates=['Time'], index_col='Time')
df=df[(df.index >= '2024-10-08 23:00:00') & (df.index < '2024-10-09 02:30:00')]  # Filter data after 00:00 UTC on 2024-11-01
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
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.tight_layout()
plt.savefig('goes_xray_lightcurve.png', dpi=300)
plt.show()


