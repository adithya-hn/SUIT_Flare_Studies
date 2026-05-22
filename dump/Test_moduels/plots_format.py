import numpy as np
import matplotlib.pyplot as plt

# Set up matplotlib to use LaTeX for text rendering
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "figure.titlesize": 16,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.format": "pdf",
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.1,
    "lines.linewidth": 2,
    "axes.linewidth": 1.5,
    "xtick.major.width": 1.5,
    "ytick.major.width": 1.5,
    "xtick.minor.width": 1,
    "ytick.minor.width": 1,
    "xtick.major.size": 6,
    "ytick.major.size": 6,
    "xtick.minor.size": 3,
    "ytick.minor.size": 3,
    "xtick.direction": "in",  # Ticks point inward
    "ytick.direction": "in",  # Ticks point inward
})

# Generate some example data
x = np.linspace(0, 10, 100)
y1 = np.sin(x)
y2 = np.cos(x)
y3 = np.sin(x) * np.cos(x)

# Create the plot
fig, ax = plt.subplots(figsize=(8, 6))

# Plot the data
ax.plot(x, y1, label=r"$\sin(x)$", color="blue", linestyle="-")
ax.plot(x, y2, label=r"$\cos(x)$", color="red", linestyle="--")
ax.plot(x, y3, label=r"$\sin(x) \cdot \cos(x)$", color="green", linestyle="-.")

# Customize the plot
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.set_title("Example Publication-Quality Plot")
ax.legend(loc="upper right", frameon=False)

# Remove grid
ax.grid(False)

# Ensure ticks point inward
ax.tick_params(axis='both', direction='in', which='both')

# Add minor ticks
ax.minorticks_on()

# Save the plot
plt.savefig("publication_quality_plot.pdf")
plt.show()