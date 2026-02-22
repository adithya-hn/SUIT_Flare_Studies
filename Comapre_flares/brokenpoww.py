



import numpy as np
from scipy.optimize import minimize

import matplotlib.pyplot as plt


def neg_log_likelihood(params, E, Emin):
    alpha1, alpha2, Eb = params

    if alpha1 <= 1 or alpha2 <= 1 or Eb <= Emin:
        return np.inf

    # Split data
    E1 = E[(E >= Emin) & (E < Eb)]
    E2 = E[E >= Eb]

    # Normalization
    Cinv = (
        (Eb**(1-alpha1) - Emin**(1-alpha1)) / (1-alpha1)
        + Eb**(alpha2-alpha1) * (Eb**(1-alpha2)) / (alpha2-1)
    )
    C = 1 / Cinv

    logL = 0.0
    if len(E1) > 0:
        logL += np.sum(np.log(C) - alpha1*np.log(E1))
    if len(E2) > 0:
        logL += np.sum(np.log(C) + (alpha2-alpha1)*np.log(Eb)
                       - alpha2*np.log(E2))

    return -logL


def broken_ccdf(E, Emin, alpha1, alpha2, Eb, N0):
    ccdf = np.zeros_like(E, dtype=float)

    mask1 = (E >= Emin) & (E < Eb)
    mask2 = E >= Eb

    ccdf[mask1] = N0 * (E[mask1]/Emin)**(-(alpha1-1))
    ccdf[mask2] = (
        N0 * (Eb/Emin)**(-(alpha1-1))
        * (E[mask2]/Eb)**(-(alpha2-1))
    )

    return ccdf

Efit = E[E >= Emin_opt]

init = [2.0, 2.5, np.median(Efit)]
bounds = [(1.01, 5), (1.01, 5), (Emin_opt*1.1, Efit.max())]

res = minimize(neg_log_likelihood,
               init,
               args=(Efit, Emin_opt),
               bounds=bounds)

alpha1_mle, alpha2_mle, Eb_mle = res.x

print(alpha1_mle, alpha2_mle, Eb_mle)

E_sorted = np.sort(E)
N_emp = np.arange(len(E_sorted), 0, -1)

N0 = np.sum(E >= Emin_opt)
E_plot = E_sorted[E_sorted >= Emin_opt]

ccdf_model = broken_ccdf(E_plot, Emin_opt,
                          alpha1_mle, alpha2_mle, Eb_mle, N0)

plt.figure(figsize=(6,5))
plt.loglog(E_sorted, N_emp, 'o', ms=4, label='All events')
plt.loglog(E_plot, ccdf_model, 'r-', lw=2,
           label='Broken power-law MLE')

plt.axvline(Emin_opt, ls='--', c='k', label=r'$E_{\min}$')
plt.axvline(Eb_mle, ls=':', c='r', label=r'$E_b$')

plt.xlabel('Peak excess Mg II h counts')
plt.ylabel('N(>E)')
plt.legend()
plt.grid(True, which='both', alpha=0.3)
plt.show()
