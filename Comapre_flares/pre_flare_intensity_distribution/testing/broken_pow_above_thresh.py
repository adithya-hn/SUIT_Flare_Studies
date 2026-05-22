"""
@Author      : Adithya H N
@Created On  : 2026-04-09
@Last Updated: 2026-04-09
@Project     : Pre-flare transienst
@Version     : 1.0

@Description
-----------
Brief description: Fitting MLE based broken powerlaw
"""

import numpy as np



data = np.loadtxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Comapre_flares/pre_flare_intensity_distribution/Flux_distrib/Peak_Flux.csv')
E = np.sort(data)
E=E[E>0.4]
N = len(E)


def trunc_norm(alpha, xmin, xmax):
    if np.isclose(alpha, 1.0):
        return np.log(xmax / xmin)
    return (xmax**(1 - alpha) - xmin**(1 - alpha)) / (1 - alpha)

def neg_logL_threePL(params, data, xmin, xmax):
    alpha1, alpha2, alpha3, xb1, xb2 = params

    # constraints
    if (alpha1 <= 1) or (alpha2 <= 1) or (alpha3 <= 1):
        return np.inf
    if not (xmin < xb1 < xb2 < xmax):
        return np.inf

    # split data
    seg1 = data[(data >= xmin) & (data < xb1)]
    seg2 = data[(data >= xb1) & (data < xb2)]
    seg3 = data[(data >= xb2) & (data <= xmax)]

    if len(seg1) == 0 or len(seg2) == 0 or len(seg3) == 0:
        return np.inf

    N1, N2, N3 = len(seg1), len(seg2), len(seg3)

    # normalization
    C1 = trunc_norm(alpha1, xmin, xb1)
    C2 = trunc_norm(alpha2, xb1, xb2)
    C3 = trunc_norm(alpha3, xb2, xmax)

    # negative log-likelihood
    nll = (
        alpha1 * np.sum(np.log(seg1)) + N1 * np.log(C1) +
        alpha2 * np.sum(np.log(seg2)) + N2 * np.log(C2) +
        alpha3 * np.sum(np.log(seg3)) + N3 * np.log(C3)
    )

    return nll

from scipy.optimize import minimize

def fit_three_powerlaw(data, xmin, xmax):

    # initial guesses
    init_alpha1 = 1.51
    init_alpha2 = 2.1
    init_alpha3 = 3.01

    # log-space breaks (important!)
    log_data = np.log10(data)
    xb1 = 10**np.percentile(log_data, 30)
    xb2 = 10**np.percentile(log_data, 66)

    x0 = [init_alpha1, init_alpha2, init_alpha3, xb1, xb2]

    bounds = [
        (1.01, 2),   # alpha1
        (2, 3),   # alpha2
        (2, 4),   # alpha3
        (xmin * 1.01, xmax * 0.8),  # xb1
        (xmin * 1.2, xmax * 0.99)   # xb2
    ]

    res = minimize(
        neg_logL_threePL,
        x0=x0,
        args=(data, xmin, xmax),
        bounds=bounds,
        method='L-BFGS-B'
    )

    return res.x, res
eb1=0
xmin = eb1
xmax = max(E)

params, result = fit_three_powerlaw(E, xmin, xmax)

alpha1, alpha2, alpha3, xb1, xb2 = params

print("alpha1 =", alpha1)
print("alpha2 =", alpha2)
print("alpha3 =", alpha3)
print("xb1 =", xb1)
print("xb2 =", xb2)


'''
def trunc_norm(alpha, xmin, xmax):
    if np.isclose(alpha, 1.0):
        return np.log(xmax / xmin)
    return (xmax**(1 - alpha) - xmin**(1 - alpha)) / (1 - alpha)

def neg_logL_broken(params, data, xmin, xmax):
    alpha1, alpha2, xb = params

    # constraints
    if (alpha1 <= 1) or (alpha2 <= 1):
        return np.inf
    if not (xmin < xb < xmax):
        return np.inf

    # split data
    seg1 = data[(data >= xmin) & (data < xb)]
    seg2 = data[(data >= xb) & (data <= xmax)]

    if len(seg1) == 0 or len(seg2) == 0:
        return np.inf

    N1 = len(seg1)
    N2 = len(seg2)

    # normalization constants
    C1 = trunc_norm(alpha1, xmin, xb)
    C2 = trunc_norm(alpha2, xb, xmax)

    # negative log-likelihood
    nll = (alpha1 * np.sum(np.log(seg1)) + N1 * np.log(C1) +
        alpha2 * np.sum(np.log(seg2)) + N2 * np.log(C2)
    )

    return nll

from scipy.optimize import minimize

def fit_broken_powerlaw(data, xmin, xmax):
    
    # initial guesses
    init_alpha1 = 1.5
    init_alpha2 = 2.5
    init_xb = np.median(data)
    bounds = [
        (1.01, 10),     # alpha1
        (1.01, 10),     # alpha2
        (xmin * 1.01, xmax * 0.99)  # break
    ]

    res = minimize(
        neg_logL_broken,
        x0=[init_alpha1, init_alpha2, init_xb],
        args=(data, xmin, xmax),
        bounds=bounds,
        method='L-BFGS-B'
    )

    return res.x, res

eb1=0.99
xmin = eb1
xmax = max(E)

(params, result) = fit_broken_powerlaw(E, xmin, xmax)

alpha1, alpha2, xb = params

print("alpha1 =", alpha1)
print("alpha2 =", alpha2)
print("break xb =", xb)'''