import numpy as np
from scipy.optimize import minimize

data = np.loadtxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Comapre_flares/pre_flare_intensity_distribution/Flux_distrib/Peak_Flux.csv')

data=np.sort(data)

def neg_logL_broken(params, data, xmin, xmax):
    alpha1, alpha2, xb = params

    # constraints
    if alpha1 <= 1 or alpha2 <= 1:
        return np.inf
    if xb <= xmin or xb >= xmax:
        return np.inf

    # split data
    data1 = data[data < xb]
    data2 = data[data >= xb]

    N = len(data)

    # normalization constants
    C1 = (xb**(1-alpha1) - xmin**(1-alpha1)) / (1-alpha1)
    C2 = (xmax**(1-alpha2) - xb**(1-alpha2)) / (1-alpha2)

    Z = C1 + C2

    # likelihood
    term1 = alpha1 * np.sum(np.log(data1)) if len(data1) > 0 else 0
    term2 = alpha2 * np.sum(np.log(data2)) if len(data2) > 0 else 0

    return term1 + term2 + N * np.log(Z)


def fit_broken_powerlaw(data, xmin, xmax):
    initial_guess = [2.0, 2.5, np.median(data)]

    bounds = [(1.01, 10), (1.01, 10), (xmin*1.1, xmax*0.9)]

    res = minimize(
        neg_logL_broken,
        initial_guess,
        args=(data, xmin, xmax),
        bounds=bounds,
        method='L-BFGS-B'
    )

    return res.x  # alpha1, alpha2, xb

print(fit_broken_powerlaw(data,0.36,max(data)))