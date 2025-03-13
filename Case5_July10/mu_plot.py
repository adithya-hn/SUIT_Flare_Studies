import numpy as np
import matplotlib.pyplot as plt

# Solar parameters (in arcseconds)
R_sun_arcsec = 960  # Approximate solar radius in arcseconds
sun_center_x, sun_center_y = 0, 0  # Sun center at (0,0) in HPC

# Example flare locations (HPC in arcsec) and their flare classes
flare_data = [
     (-438, -291,"X1.4"),  # (X, Y) in arcseconds
    (-263, 298,"M1.2"),
    (-210, -285,"M2.1"),
    (-148,197,"M1.5"),
    (-65,-201,"M1.0"),
    (109,119,"X1.8"),
    (-413,225,"M1.3"),
    (-324,234,"M2.0"),
    (-43,-260,"M1.0"),
    (152,-206,"M1.7"),
    (264,200,"M7.6"),
]

# Compute mu values
flare_mu = []
for x_arcsec, y_arcsec, _ in flare_data:
    r = np.sqrt((x_arcsec - sun_center_x)**2 + (y_arcsec - sun_center_y)**2) / R_sun_arcsec
    mu = np.sqrt(1 - r**2) if r <= 1 else 0
    flare_mu.append(mu)

# Plotting
fig, ax = plt.subplots(figsize=(6, 6))

# Solar disk
circle = plt.Circle((sun_center_x, sun_center_y), R_sun_arcsec, color='orange', fill=False, linewidth=2)
ax.add_patch(circle)

# Scatter plot of flares
sc = ax.scatter(
    [x for x, y, _ in flare_data], 
    [y for x, y, _ in flare_data], 
    c=flare_mu, cmap='rainbow', edgecolors='black', s=50, label="Flare Locations"
)

# Annotations for each flare point
for (x, y, flare_class), mu in zip(flare_data, flare_mu):
    ax.annotate(f"{flare_class}\nμ={mu:.3f}", (x, y), textcoords="offset points", xytext=(5, 5), ha='left', fontsize=10, color='black')

# Colorbar for mu values
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label(r'$\mu$ (cos θ)')

# Labels and limits
ax.set_xlim(-R_sun_arcsec, R_sun_arcsec)
ax.set_ylim(-R_sun_arcsec, R_sun_arcsec)
ax.set_xlabel("Helioprojective X (arcsec)")
ax.set_ylabel("Helioprojective Y (arcsec)")
ax.set_title("Flare Locations on the Solar Disk")


# Equal aspect ratio
ax.set_aspect('equal')

# Show the plot
plt.show()
