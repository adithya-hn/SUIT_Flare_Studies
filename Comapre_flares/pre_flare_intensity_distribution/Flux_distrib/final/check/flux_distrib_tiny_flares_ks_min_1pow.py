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

scol =sns.color_palette("colorblind")
E=np.loadtxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Comapre_flares/pre_flare_intensity_distribution/Flux_distrib/final/Peak_Flux.csv')
E = E[E > 0]
E = np.sort(E)



def mle_alpha(E, Emin):
    return 1.0 + len(E) / np.sum(np.log(E / Emin))

def empirical_ccdf(E):
    n = len(E)
    return np.arange(n, 0, -1) / n

def model_ccdf(E, Emin, alpha):
    return (E / Emin) ** (-(alpha - 1))

def ks_distance(E, Emin):
    Efit = E[E >= Emin]
    n = len(Efit)

    # Require enough statistics
    if n < 30:
        return np.nan, np.nan

    alpha = mle_alpha(Efit, Emin)

    ccdf_emp = empirical_ccdf(Efit)
    ccdf_mod = model_ccdf(Efit, Emin, alpha)

    D = np.max(np.abs(ccdf_emp - ccdf_mod))
    return D, alpha
Emins = np.unique(E)
ks_vals = []
alphas = []

for Emin in Emins:
    D, a = ks_distance(E, Emin)
    ks_vals.append(D)
    alphas.append(a)

ks_vals = np.array(ks_vals)
alphas = np.array(alphas)
valid = ~np.isnan(ks_vals)

Emins_valid = Emins[valid]
ks_valid = ks_vals[valid]
alphas_valid = alphas[valid]

idx = np.argmin(ks_valid)

Emin_opt = Emins_valid[idx]

alpha_opt = alphas_valid[idx]

print("Optimal Emin =", Emin_opt)
print("MLE alpha =", alpha_opt)

Efit = E[E >= Emin_opt]
alpha_final = mle_alpha(Efit, Emin_opt)

print("Final alpha =", alpha_final)
print("Number of events =", len(Efit))

sigma_alpha = (alpha_final - 1) / np.sqrt(len(Efit))
print("Alpha uncertainty ~", sigma_alpha)

def get_ccdf_y(E_sorted, ccdf_full, values):
    y_vals = []
    
    for v in values:
        idx = np.searchsorted(E_sorted, v, side='left')
        if idx < len(ccdf_full):
            y_vals.append(ccdf_full[idx])
        else:
            y_vals.append(np.nan)
    
    return np.array(y_vals)

Efit = np.sort(Efit)
ccdf_emp = empirical_ccdf(Efit)

ccdf_model = model_ccdf(Efit, Emin_opt, alpha_final)
E_all = np.sort(E)                     # all base-subtracted peak counts
N_all = np.arange(len(E_all), 0, -1)
N0 = len(Efit)  # number of events above Emin
model_ccdf = 0.33 * (Efit / Emin_opt) ** (-(alpha_final - 1))

ccdf = empirical_ccdf(E)
# y_matched = get_ccdf_y(E, ccdf, matched_ints)

plt.figure(figsize=(12,8))
plt.title('CCDF of all pre-flare transients',fontsize=24)
# Full CCDF (all events)
# plt.loglog(E_all, N_all, '+',color='k', markersize=8,label='All detected events')

# Power-law fit (only above Emin)
plt.loglog(E, ccdf, '+',color='k', label='SUIT transeint events')
plt.loglog(Efit, model_ccdf, '-',
           linewidth=2,
           label=fr'Power-law fit (alpha = {alpha_final:.2f}$\pm$ {sigma_alpha:.2f} )')

# plt.scatter(matched_ints, y_matched,
#             facecolors='none', edgecolors='red',
#             s=70, linewidth=1.5,
#             label='Co-temporal HEL1OS')

# Emin marker
plt.axvline(Emin_opt, color='b', linestyle='--',label=r'$E_{\min}$')

# y_matched = get_ccdf_y(E, ccdf, matched_ints)

plt.title('CCDF of pre-flare transients')
plt.xlabel(r'Energy ($erg^-1 \times cm^-2$)', fontsize=18)
plt.ylabel('CCDF', fontsize=18)
plt.legend(fontsize=14)
plt.grid(True, which='both', alpha=0.3)
plt.savefig('ccdf_mle_ks_min.pdf',dpi=300)
plt.show()

def ks_pvalue(E, Emin, alpha_obs, n_synthetic=1000):
    """
    Clauset et al. 2009 Section 4 — synthetic dataset p-value.
    p > 0.1 means power law is a plausible fit.
    """
    Efit = E[E >= Emin]
    n    = len(Efit)

    # Observed KS distance
    ccdf_emp = np.arange(n, 0, -1) / n
    ccdf_mod = (Efit / Emin) ** (-(alpha_obs - 1))
    D_obs    = np.max(np.abs(ccdf_emp - ccdf_mod))

    # Synthetic datasets drawn from the fitted power law
    rng     = np.random.default_rng(42)
    D_synth = []
    for _ in range(n_synthetic):
        # Generate synthetic power-law data
        u       = rng.uniform(size=n)
        x_synth = Emin * (1 - u) ** (-1 / (alpha_obs - 1))
        x_synth = np.sort(x_synth)

        # Refit alpha on synthetic data
        a_s     = mle_alpha(x_synth, Emin)
        mod_s   = (x_synth / Emin) ** (-(a_s - 1))
        emp_s   = np.arange(n, 0, -1) / n
        D_synth.append(np.max(np.abs(emp_s - mod_s)))

    p = np.mean(np.array(D_synth) >= D_obs)
    return p, D_obs

p_val, D_obs = ks_pvalue(E, Emin_opt, alpha_final, n_synthetic=2000)
print(f"KS D = {D_obs:.4f}")
print(f"p-value = {p_val:.3f}")
print("Power law plausible" if p_val > 0.1 else "Power law rejected")

'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.stats import kstest

# ── Your data ─────────────────────────────────────────────────────────────────
data = E
N    = len(data)
ccdf = (N - np.arange(N)) / N

# ── Set xmin ABOVE the incompleteness plateau ─────────────────────────────────
# This is the most important choice — set it where the CCDF first
# begins falling consistently (visually ~0.07 from your plot)
xmin = 0.36        # tune this — everything below is discarded

data_fit = data[data >= xmin]
ccdf_fit = ccdf[data >= xmin]     # matching empirical CCDF values

# ── Single breakpoint between the two real power laws ─────────────────────────
b = 1.0           # from your visual inspection

seg2 = data_fit[data_fit <  b]
seg3 = data_fit[data_fit >= b]

print(f"N above xmin : {len(data_fit)}")
print(f"N in Seg 2   : {len(seg2)}")
print(f"N in Seg 3   : {len(seg3)}")

# ── MLE Hill estimator ────────────────────────────────────────────────────────
def mle_alpha(x, xmin_seg):
    x = x[x >= xmin_seg]
    n = len(x)
    if n < 5:
        return np.nan, -np.inf
    alpha = 1.0 + n / np.sum(np.log(x / xmin_seg))
    logL  = (n * np.log(alpha - 1) - n * np.log(xmin_seg)
             - alpha * np.sum(np.log(x / xmin_seg)))
    return alpha, logL

a2, L2 = mle_alpha(seg2, xmin)
a3, L3 = mle_alpha(seg3, b)
print(f"\nMLE:  α2 = {a2:.3f}   α3 = {a3:.3f}")

# ── Bootstrap uncertainties ───────────────────────────────────────────────────
def bootstrap_alpha(seg, xmin_seg, n_boot=2000):
    rng    = np.random.default_rng(42)
    alphas = []
    for _ in range(n_boot):
        s = np.sort(rng.choice(seg, len(seg), replace=True))
        a, _ = mle_alpha(s, xmin_seg)
        if not np.isnan(a):
            alphas.append(a)
    return np.mean(alphas), np.std(alphas)

a2m, a2s = bootstrap_alpha(seg2, xmin)
a3m, a3s = bootstrap_alpha(seg3, b)
print(f"Bootstrap: α2 = {a2m:.3f} ± {a2s:.3f}")
print(f"           α3 = {a3m:.3f} ± {a3s:.3f}")

# ── KS test ───────────────────────────────────────────────────────────────────
def ks_pl(seg, alpha, xmin_seg):
    s   = seg[seg >= xmin_seg]
    cdf = lambda x: 1 - (xmin_seg / x) ** (alpha - 1)
    D, p = kstest(s, cdf)
    return D, p

D2, p2 = ks_pl(seg2, a2m, xmin)
D3, p3 = ks_pl(seg3, a3m, b)
print(f"\nKS test:  seg2 D={D2:.3f} p={p2:.3f}")
print(f"          seg3 D={D3:.3f} p={p3:.3f}")

# ── Anchored fit lines (log-space least squares on CCDF) ──────────────────────
def anchor(x_seg, y_seg, alpha):
    log_C = np.mean(np.log(y_seg) - (1 - alpha) * np.log(x_seg))
    return np.exp(log_C)

m2 = (data_fit >= xmin) & (data_fit <  b)
m3 =  data_fit >= b

C2 = anchor(data_fit[m2], ccdf_fit[m2], a2m)
C3 = anchor(data_fit[m3], ccdf_fit[m3], a3m)

def pl(x, alpha, C):
    return C * x ** (1 - alpha)

x2_fit = np.logspace(np.log10(xmin), np.log10(b),           300)
x3_fit = np.logspace(np.log10(b),    np.log10(data[-1]*1.05), 300)

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 1, figsize=(9, 9),
                          gridspec_kw={"height_ratios": [3, 1]},
                          sharex=True)
ax, ax_res = axes

# All data (including incomplete region) shown as faint crosses
ax.scatter(data[data < xmin], ccdf[data < xmin],
           s=12, marker="+", color="gray", alpha=0.4,
           label="Below completeness limit", zorder=2)

# Fitted data
ax.scatter(data_fit, ccdf_fit,
           s=15, marker="+", color="black",
           label="SUIT transient events", zorder=3)

ax.plot(x2_fit, pl(x2_fit, a2m, C2), color="#E07B39", lw=2.2,
        label=fr"$\alpha_1$ = {a2m:.2f} ± {a2s:.2f}  (KS p={p2:.2f})")
ax.plot(x3_fit, pl(x3_fit, a3m, C3), color="#2D9E55", lw=2.2,
        label=fr"$\alpha_2$ = {a3m:.2f} ± {a3s:.2f}  (KS p={p3:.2f})")

# Completeness limit
ax.axvline(xmin, color="gray",    lw=1.0, ls=":",  alpha=0.8,
           label=f"Completeness limit ({xmin})")
# Breakpoint
ax.axvline(b,    color="#E07B39", lw=0.8, ls="--", alpha=0.6,
           label=f"Break energy ({b} erg cm⁻²)")

ax.set_yscale("log"); ax.set_xscale("log")
ax.set_ylabel("CCDF", fontsize=12)
ax.legend(fontsize=9, framealpha=0.9)
ax.set_title("CCDF of pre-flare transients — MLE power-law fit", fontsize=13)
ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)

# ── Residuals (fitted region only) ───────────────────────────────────────────
def model(x):
    return np.where(x < b, pl(x, a2m, C2), pl(x, a3m, C3))

res = (ccdf_fit - model(data_fit)) / model(data_fit)

ax_res.scatter(data_fit, res, s=10, marker="+", color="black")
ax_res.axhline(0,    color="gray", lw=0.8)
ax_res.axhline( 0.2, color="gray", lw=0.5, ls="--")
ax_res.axhline(-0.2, color="gray", lw=0.5, ls="--")
ax_res.axvline(b,    color="gray", lw=0.8, ls="--", alpha=0.5)

ax_res.set_ylabel("Frac. residual", fontsize=11)
ax_res.set_xlabel(r"Energy (erg$^{-1}$ × cm$^{-2}$)", fontsize=12)
ax_res.set_ylim(-0.6, 0.6)
ax_res.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)

plt.tight_layout()
plt.savefig("ccdf_mle_final.pdf", dpi=180, bbox_inches="tight")
plt.show()'''

# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.ticker as ticker

# # ── 0. Your data ──────────────────────────────────────────────────────────────
# data =E# np.sort(data)
# N    = len(data)
# ccdf = (N - np.arange(N)) / N

# # ── 1. SET BREAKPOINTS MANUALLY from visual inspection ────────────────────────
# # Tweak these two values until the lines sit on the data
# b1 = 0.36    # end of flat region / start of mid slope
# b2 = 1.0    # start of steep drop

# seg1 = data[data < b1]
# seg2 = data[(data >= b1) & (data < b2)]
# seg3 = data[data >= b2]

# print(f"N per segment: {len(seg1)}, {len(seg2)}, {len(seg3)}")

# # ── 2. MLE α per segment ──────────────────────────────────────────────────────
# def mle_alpha(x, xmin):
#     x = x[x >= xmin]
#     n = len(x)
#     if n < 4:
#         return np.nan, -np.inf
#     alpha = 1.0 + n / np.sum(np.log(x / xmin))
#     logL  = n * np.log(alpha - 1) - n * np.log(xmin) \
#             - alpha * np.sum(np.log(x / xmin))
#     return alpha, logL

# a1, _ = mle_alpha(seg1, seg1[0])
# a2, _ = mle_alpha(seg2, b1)
# a3, _ = mle_alpha(seg3, b2)
# print(f"α1={a1:.3f}  α2={a2:.3f}  α3={a3:.3f}")

# # ── 3. Bootstrap uncertainty ──────────────────────────────────────────────────
# def bootstrap_alpha(seg, xmin, n_boot=1000):
#     rng = np.random.default_rng(42)
#     alphas = [mle_alpha(np.sort(rng.choice(seg, len(seg), replace=True)), xmin)[0]
#               for _ in range(n_boot)]
#     alphas = [a for a in alphas if not np.isnan(a)]
#     return np.mean(alphas), np.std(alphas)

# a1m, a1s = bootstrap_alpha(seg1, seg1[0])
# a2m, a2s = bootstrap_alpha(seg2, b1)
# a3m, a3s = bootstrap_alpha(seg3, b2)
# print(f"Bootstrap: α1={a1m:.3f}±{a1s:.3f}  α2={a2m:.3f}±{a2s:.3f}  α3={a3m:.3f}±{a3s:.3f}")

# # ── 4. Anchored power-law lines ───────────────────────────────────────────────
# def anchor(data_full, ccdf_full, alpha, x_lo, x_hi):
#     mask  = (data_full >= x_lo) & (data_full < x_hi)
#     x_seg = data_full[mask]
#     y_seg = ccdf_full[mask]
#     if len(x_seg) < 2:
#         return 1.0
#     log_C = np.mean(np.log(y_seg) - (1 - alpha) * np.log(x_seg))
#     return np.exp(log_C)

# C1 = anchor(data, ccdf, a1m, seg1[0], b1)
# C2 = anchor(data, ccdf, a2m, b1,      b2)
# C3 = anchor(data, ccdf, a3m, b2,      data[-1])

# def pl(x, alpha, C):
#     return C * x ** (1 - alpha)

# x1 = np.logspace(np.log10(max(data[0],  1e-3)), np.log10(b1),        300)
# x2 = np.logspace(np.log10(b1),                  np.log10(b2),        300)
# x3 = np.logspace(np.log10(b2),                  np.log10(data[-1]),  300)

# # ── 5. Plot ───────────────────────────────────────────────────────────────────
# fig, axes = plt.subplots(2, 1, figsize=(9, 9),
#                           gridspec_kw={"height_ratios": [3, 1]},
#                           sharex=True)
# ax, ax_res = axes

# ax.scatter(data, ccdf, s=15, marker="+", color="black",
#            label="SUIT transient events", zorder=3)

# ax.plot(x1, pl(x1, a1m, C1), color="#2E86AB", lw=2,
#         label=f"Seg 1  α={a1m:.2f}±{a1s:.2f}")
# ax.plot(x2, pl(x2, a2m, C2), color="#E07B39", lw=2,
#         label=f"Seg 2  α={a2m:.2f}±{a2s:.2f}")
# ax.plot(x3, pl(x3, a3m, C3), color="#2D9E55", lw=2,
#         label=f"Seg 3  α={a3m:.2f}±{a3s:.2f}")

# for bk, c in [(b1, "#2E86AB"), (b2, "#E07B39")]:
#     ax.axvline(bk, color=c, lw=0.8, ls="--", alpha=0.6)

# ax.set_yscale("log"); ax.set_xscale("log")
# ax.set_ylabel("CCDF", fontsize=12)
# ax.legend(fontsize=10, framealpha=0.9)
# ax.set_title("CCDF of pre-flare transients — piecewise power-law MLE fit", fontsize=13)
# ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)

# # ── 6. Residuals ──────────────────────────────────────────────────────────────
# def model(x):
#     return np.where(x < b1, pl(x, a1m, C1),
#            np.where(x < b2, pl(x, a2m, C2),
#                              pl(x, a3m, C3)))

# res = (ccdf - model(data)) / model(data)
# ax_res.scatter(data, res, s=10, marker="+", color="black")
# ax_res.axhline(0,    color="gray", lw=0.8)
# ax_res.axhline( 0.2, color="gray", lw=0.5, ls="--")
# ax_res.axhline(-0.2, color="gray", lw=0.5, ls="--")
# for bk in [b1, b2]:
#     ax_res.axvline(bk, color="gray", lw=0.8, ls="--", alpha=0.5)
# ax_res.set_ylabel("Frac. residual", fontsize=11)
# ax_res.set_xlabel(r"Energy (erg$^{-1}$ × cm$^{-2}$)", fontsize=12)
# ax_res.set_ylim(-0.6, 0.6)
# ax_res.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)

# plt.tight_layout()
# plt.savefig("ccdf_mle_fit_v2.pdf", dpi=180, bbox_inches="tight")
# plt.show()

#------------------------------------------------------------------------------------------------------------------------
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.ticker as ticker
# from scipy.stats import kstest
# from itertools import product

# # ── 0. Load your data ────────────────────────────────────────────────────────
# # Replace this with your actual energy array
# # data = np.load("your_energies.npy")          # shape (N,), units: erg cm⁻²
# data = E

# # ── 1. MLE (Hill estimator) for a single power-law segment ───────────────────
# def mle_alpha(x, xmin):
#     x = x[x >= xmin]
#     n = len(x)
#     if n < 6:
#         return np.nan, -np.inf
#     alpha = 1.0 + n / np.sum(np.log(x / xmin))
#     logL  = n * np.log(alpha - 1) - n * np.log(xmin) \
#             - alpha * np.sum(np.log(x / xmin))
#     return alpha, logL

# # ── 2. Grid-search for the two best breakpoints ───────────────────────────────
# def find_breakpoints(data, n_grid=60, min_seg=8):
#     lo, hi = data[min_seg], data[-(min_seg+1)]
#     cands  = np.logspace(np.log10(lo), np.log10(hi), n_grid)
#     best_L, best = -np.inf, None

#     for b1, b2 in product(cands, cands):
#         if b1 >= b2:
#             continue
#         s1 = data[data <  b1]
#         s2 = data[(data >= b1) & (data < b2)]
#         s3 = data[data >= b2]
#         if min(len(s1), len(s2), len(s3)) < min_seg:
#             continue
#         a1, L1 = mle_alpha(s1, s1[0])
#         a2, L2 = mle_alpha(s2, b1)
#         a3, L3 = mle_alpha(s3, b2)
#         if np.isnan(a1+a2+a3):
#             continue
#         total = L1 + L2 + L3
#         if total > best_L:
#             best_L, best = total, (b1, b2, a1, a2, a3)

#     return best   # (b1, b2, α1, α2, α3)

# b1, b2, a1, a2, a3 = find_breakpoints(data)
# print(f"Breakpoints : {b1:.4f}, {b2:.4f}")
# print(f"α₁ = {a1:.3f}  |  α₂ = {a2:.3f}  |  α₃ = {a3:.3f}")

# # ── 3. Bootstrap uncertainty on each α ───────────────────────────────────────
# def bootstrap_alpha(seg, xmin, n_boot=1000, rng=None):
#     rng = rng or np.random.default_rng(42)
#     alphas = []
#     for _ in range(n_boot):
#         s      = rng.choice(seg, size=len(seg), replace=True)
#         a, _   = mle_alpha(np.sort(s), xmin)
#         if not np.isnan(a):
#             alphas.append(a)
#     return np.mean(alphas), np.std(alphas)

# seg1 = data[data < b1]
# seg2 = data[(data >= b1) & (data < b2)]
# seg3 = data[data >= b2]

# rng = np.random.default_rng(42)
# a1m, a1s = bootstrap_alpha(seg1, seg1[0],  rng=rng)
# a2m, a2s = bootstrap_alpha(seg2, b1,        rng=rng)
# a3m, a3s = bootstrap_alpha(seg3, b2,        rng=rng)
# print(f"\nBootstrap:")
# print(f"  α₁ = {a1m:.3f} ± {a1s:.3f}")
# print(f"  α₂ = {a2m:.3f} ± {a2s:.3f}")
# print(f"  α₃ = {a3m:.3f} ± {a3s:.3f}")

# # ── 4. KS test per segment ────────────────────────────────────────────────────
# def ks_test_powerlaw(seg, alpha, xmin):
#     cdf = lambda x: 1 - (xmin / x) ** (alpha - 1)
#     D, p = kstest(seg[seg >= xmin], cdf)
#     return D, p

# D1, p1 = ks_test_powerlaw(seg1, a1, seg1[0])
# D2, p2 = ks_test_powerlaw(seg2, a2, b1)
# D3, p3 = ks_test_powerlaw(seg3, a3, b2)
# print(f"\nKS test  (p > 0.1 = good fit):")
# print(f"  seg1: D={D1:.3f}, p={p1:.3f}")
# print(f"  seg2: D={D2:.3f}, p={p2:.3f}")
# print(f"  seg3: D={D3:.3f}, p={p3:.3f}")

# # ── 5. Helper: theoretical CCDF line for one power-law segment ─────────────────
# def powerlaw_ccdf(x_range, alpha, xmin, norm):
#     """norm = CCDF value at xmin (so fits attach to the data)."""
#     return norm * (xmin / x_range) ** (alpha - 1)

# # Empirical CCDF
# N    = len(data)
# ccdf = (N - np.arange(N)) / N      # P(X > x_i)

# # Anchor normalisation at each breakpoint
# norm1 = np.interp(seg1[0], data, ccdf)
# norm2 = np.interp(b1,      data, ccdf)
# norm3 = np.interp(b2,      data, ccdf)

# x1_fit = np.logspace(np.log10(seg1[0]),  np.log10(b1),  200)
# x2_fit = np.logspace(np.log10(b1),       np.log10(b2),  200)
# x3_fit = np.logspace(np.log10(b2),       np.log10(data[-1]*1.05), 200)

# # ── 6. Main plot: CCDF + fits ─────────────────────────────────────────────────
# fig, axes = plt.subplots(2, 1, figsize=(9, 9),
#                           gridspec_kw={"height_ratios": [3, 1]},
#                           sharex=True)
# ax, ax_res = axes

# # Data
# ax.scatter(data, ccdf, s=15, marker="+", color="black",
#            label="SUIT transient events", zorder=3)

# # Fitted power-law lines
# kw = dict(lw=2.0, alpha=0.85)
# ax.plot(x1_fit, powerlaw_ccdf(x1_fit, a1m, seg1[0], norm1),
#         color="#2E86AB", label=f"Seg 1  α={a1m:.2f}±{a1s:.2f}", **kw)
# ax.plot(x2_fit, powerlaw_ccdf(x2_fit, a2m, b1,       norm2),
#         color="#E07B39", label=f"Seg 2  α={a2m:.2f}±{a2s:.2f}", **kw)
# ax.plot(x3_fit, powerlaw_ccdf(x3_fit, a3m, b2,       norm3),
#         color="#2D9E55", label=f"Seg 3  α={a3m:.2f}±{a3s:.2f}", **kw)

# # Breakpoint verticals
# for bk, col in [(b1, "#2E86AB"), (b2, "#2D9E55")]:
#     ax.axvline(bk, color=col, lw=0.8, ls="--", alpha=0.6)

# ax.set_yscale("log");  ax.set_xscale("log")
# ax.set_ylabel("CCDF", fontsize=12)
# ax.legend(fontsize=10, framealpha=0.9)
# ax.set_title("CCDF of pre-flare transients — piecewise power-law MLE fit",
#              fontsize=13)
# ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)
# ax.xaxis.set_major_formatter(ticker.LogFormatterMathtext())
# ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext())

# # ── 7. Residual panel: (data − model) / model ────────────────────────────────
# def model_ccdf(x):
#     if x < b1:
#         return powerlaw_ccdf(np.array([x]), a1m, seg1[0], norm1)[0]
#     elif x < b2:
#         return powerlaw_ccdf(np.array([x]), a2m, b1,       norm2)[0]
#     else:
#         return powerlaw_ccdf(np.array([x]), a3m, b2,       norm3)[0]

# model_vals = np.array([model_ccdf(xi) for xi in data])
# residuals  = (ccdf - model_vals) / model_vals   # fractional residual

# ax_res.scatter(data, residuals, s=10, marker="+", color="black")
# ax_res.axhline(0, color="gray", lw=0.8)
# ax_res.axhline( 0.2, color="gray", lw=0.5, ls="--")
# ax_res.axhline(-0.2, color="gray", lw=0.5, ls="--")
# for bk in [b1, b2]:
#     ax_res.axvline(bk, color="gray", lw=0.8, ls="--", alpha=0.5)
# ax_res.set_ylabel("Frac. residual", fontsize=11)
# ax_res.set_xlabel(r"Energy (erg$^{-1}$ × cm$^{-2}$)", fontsize=12)
# ax_res.set_ylim(-0.6, 0.6)
# ax_res.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)

# plt.tight_layout()
# plt.savefig("ccdf_mle_fit.pdf", dpi=180, bbox_inches="tight")
# plt.show()