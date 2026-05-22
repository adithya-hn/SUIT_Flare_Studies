"""
@Author      : Adithya H N
@Created On  : 2026-04-09
@Last Updated: 2026-04-09
@Project     : Pre-flare study
@Version     : 1.0

@Description
-----------
Brief description: check ccdf to alpha relation
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import matplotlib.dates as mdates
import pandas as pd
from scipy.signal import find_peaks
from scipy.signal import argrelextrema
import numpy.ma as ma
import seaborn as sns
import glob
import os
from scipy.optimize import curve_fit
import scienceplots
plt.style.use('science')
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution, curve_fit
from scipy.stats import kstest
from scipy.optimize import minimize_scalar
from itertools import product
from scipy.optimize import minimize

data = np.loadtxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Comapre_flares/pre_flare_intensity_distribution/Flux_distrib/Peak_Flux.csv')
E=np.sort(data[data>0.33])

eb1 = 0.36
eb2 = 0.99

seg1 = E[(E >= eb1) & (E < eb2)]
seg2 = E[E >= eb2]

def mle_alpha(E, Emin):
    return 1.0 + len(E) / np.sum(np.log(E / Emin))


def neg_logL_truncated(alpha, data, xmin, xmax):
    if alpha <= 1:
        return np.inf

    N = len(data)

    # normalization term
    norm = (xmax**(1-alpha) - xmin**(1-alpha)) / (1-alpha)

    return (   alpha * np.sum(np.log(data)) +   N * np.log(norm)  )


def fit_alpha_truncated(data, xmin, xmax):
    res = minimize_scalar(
        neg_logL_truncated,
        bounds=(1.01, 10),
        args=(data, xmin, xmax),
        method='bounded'
    )
    return res.x


def empirical_ccdf(E):
    n = len(E)
    return np.arange(n, 0, -1) / n

def model_ccdf(E, Emin, alpha):
    return (E / Emin) ** (-(alpha - 1))

def get_ccdf_y(E_sorted, ccdf_full, values):
    y_vals = []
    
    for v in values:
        idx = np.searchsorted(E_sorted, v, side='left')
        if idx < len(ccdf_full):
            y_vals.append(ccdf_full[idx])
        else:
            y_vals.append(np.nan)
    
    return np.array(y_vals)

seg1 = E[(E >= eb1) & (E < eb2)]
seg2 = E[E >= eb2]

alpha1 = fit_alpha_truncated(seg1, eb1, eb2)   # correct
alpha2 = mle_alpha(seg2, eb2)                  # your formula is OK here

print(alpha1, alpha2)

E1fit=seg1
E2fit=seg2


E_all = np.sort(E)                     # all base-subtracted peak counts
N_all = np.arange(len(E_all), 0, -1)
N0 = len(E1fit)  # number of events above Emin
model_ccdf_seg1 = 0.64 * (E1fit / eb1) ** (-(alpha1 - 1))
model_ccdf_seg2 = 0.3 * (E2fit / eb2) ** (-(alpha2 - 1))



ccdf = empirical_ccdf(E)
# y_matched = get_ccdf_y(E, ccdf, matched_ints)

alpha_sig1 = (alpha1 - 1) / np.sqrt(len(E1fit))
alpha_sig2 = (alpha2 - 1) / np.sqrt(len(E2fit))


plt.figure(figsize=(12,8))
plt.title('CCDF of all pre-flare transients',fontsize=24)

# Power-law fit (only above Emin)
plt.loglog(E, ccdf, '+',color='k', label='SUIT transeint events')
# plt.loglog(E1fit, model_ccdf_seg1, 'o',linewidth=2,label=fr'Power-law fit (alpha = {alpha1:.2f}$\pm$ {alpha_sig1:.2f} )')
# plt.loglog(E2fit, model_ccdf_seg2, 'o',linewidth=2,label=fr'Power-law fit (alpha = {alpha2:.2f}$\pm$ {alpha_sig2:.2f} )')

# plt.scatter(matched_ints, y_matched,  facecolors='none', edgecolors='red',s=70, linewidth=1.5,  label='Co-temporal HEL1OS')

# Emin marker
plt.axvline(eb1, color='b', linestyle='--',label=r'$E_{\min}$')

# y_matched = get_ccdf_y(E, ccdf, matched_ints)

plt.title('CCDF of pre-flare transients')
plt.xlabel(r'Energy ($erg^-1 \times cm^-2$)', fontsize=18)
plt.ylabel('CCDF', fontsize=18)
plt.legend(fontsize=14)
plt.grid(True, which='both', alpha=0.3)
plt.savefig('ccdf_t.pdf',dpi=300)
plt.show()


# mod_E1=np.array(min(seg1),max,seg1,100)
# mod_E2=np.array(min(seg2),max,seg2,100)
# model_ccdf_seg1 = 0.64 * (mod_E1 / eb1) ** (-(alpha1 - 1))
# model_ccdf_seg2 = 0.3 * (mod_E2 / eb2) ** (-(alpha2 - 1))
# plt.loglog(E1fit, model_ccdf_seg1, '-',linewidth=2,label=fr'Power-law fit (alpha = {alpha1:.2f}$\pm$ {alpha_sig1:.2f} )')
# plt.loglog(E2fit, model_ccdf_seg2, '-',linewidth=2,label=fr'Power-law fit (alpha = {alpha2:.2f}$\pm$ {alpha_sig2:.2f} )')

# plt.show()
# alpha_sig1 = (alpha1 - 1) / np.sqrt(len(E1fit))
# alpha_sig2 = (alpha2 - 1) / np.sqrt(len(E2fit))
