import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# --------------------------------------------------
# INPUT DATA
# --------------------------------------------------
data = np.loadtxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Comapre_flares/pre_flare_intensity_distribution/Flux_distrib/Peak_Flux.csv')

data=np.sort(data[data>.36])
data=data[data<1.68]

xmin = data.min()
xmax = data.max()


# --------------------------------------------------
# Broken Power-Law PDF normalization
# --------------------------------------------------
def broken_norm(alpha1, alpha2, xb, xmin, xmax):

    I1 = (xb**(1-alpha1) - xmin**(1-alpha1)) / (1-alpha1)

    I2 = xb**(alpha2-alpha1) * \
         (xmax**(1-alpha2) - xb**(1-alpha2)) / (1-alpha2)

    return I1 + I2


# --------------------------------------------------
# Negative Log-Likelihood
# --------------------------------------------------
def neg_logL(params, data, xmin, xmax):

    alpha1, alpha2, xb = params

    if alpha1 <= 1 or alpha2 <= 1:
        return np.inf

    if not (xmin < xb < xmax):
        return np.inf

    norm = broken_norm(alpha1, alpha2, xb, xmin, xmax)

    seg1 = data[data < xb]
    seg2 = data[data >= xb]

    ll1 = -alpha1 * np.sum(np.log(seg1))

    ll2 = np.sum(
        np.log(xb**(alpha2-alpha1)) - alpha2*np.log(seg2)
    )

    N = len(data)

    return -(ll1 + ll2 - N*np.log(norm))


# --------------------------------------------------
# Optimization
# --------------------------------------------------
guess = [1., 3.0, np.median(data)]

res = minimize(
    neg_logL,
    guess,
    args=(data, xmin, xmax),
    method='Nelder-Mead'
)

alpha1, alpha2, xb = res.x

print("Best-fit parameters:")
print(f"alpha1 = {alpha1:.3f}")
print(f"alpha2 = {alpha2:.3f}")
print(f"Break  = {xb:.3f}")


# --------------------------------------------------
# Plot
# --------------------------------------------------
bins = np.logspace(np.log10(xmin), np.log10(xmax), 25)
hist, edges = np.histogram(data, bins=bins, density=True)
centers = np.sqrt(edges[:-1]*edges[1:])

plt.figure(figsize=(6,5))
plt.loglog(centers, hist, 'o', label='Data')

x1 = np.logspace(np.log10(xmin), np.log10(xb), 100)
x2 = np.logspace(np.log10(xb), np.log10(xmax), 100)

norm = 1 / broken_norm(alpha1, alpha2, xb, xmin, xmax)

y1 = norm * x1**(-alpha1)
y2 = norm * xb**(alpha2-alpha1) * x2**(-alpha2)

plt.loglog(x1, y1, '-', lw=2)
plt.loglog(x2, y2, '-', lw=2, label='Broken PL Fit')

plt.axvline(xb, ls='--', color='k', label='Break')

plt.xlabel("Energy")
plt.ylabel("PDF")
plt.legend()
plt.show()