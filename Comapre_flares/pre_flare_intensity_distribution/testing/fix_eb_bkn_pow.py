import numpy as np
from scipy.optimize import minimize

data = np.loadtxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Comapre_flares/pre_flare_intensity_distribution/Flux_distrib/Peak_Flux.csv')

data=np.sort(data[data>0.36])

xmin = data.min()
xmax = data.max()

xb = 1.0   # <-- fixed known break


# --------------------------------------------------
# Normalization
# --------------------------------------------------
def broken_norm(alpha1, alpha2, xb, xmin, xmax):

    I1 = (xb**(1-alpha1) - xmin**(1-alpha1)) / (1-alpha1)

    I2 = xb**(alpha2-alpha1) * \
         (xmax**(1-alpha2) - xb**(1-alpha2)) / (1-alpha2)

    return I1 + I2


# --------------------------------------------------
# Negative Log Likelihood
# --------------------------------------------------
def neg_logL(params):

    alpha1, alpha2 = params

    if alpha1 <= 1 or alpha2 <= 1:
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
# Optimize
# --------------------------------------------------
res = minimize(
    neg_logL,
    x0=[2.0, 3.0],
    method='Nelder-Mead'
)

alpha1, alpha2 = res.x

print(f"Fixed Break = {xb:.4f}")
print(f"Alpha1 = {alpha1:.4f}")
print(f"Alpha2 = {alpha2:.4f}")