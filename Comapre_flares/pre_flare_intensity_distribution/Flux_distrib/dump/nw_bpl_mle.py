"""
MLE Broken Power Law Fitting — Pre-Flare Transient Peak Flux
============================================================
Reference: Crawford et al., NIM-A 500 (2003)
           Clauset, Shalizi & Newman, SIAM Rev. 51 (2009)

KEY RELATIONSHIP (common source of confusion)
----------------------------------------------
  PDF  f(x) ∝ x^{-α}          → log-log slope = -α
  CCDF S(x) ∝ x^{-(α-1)}      → log-log slope = -(α-1)   ← shallower by 1!

So if you READ slope s off a log-log CCDF plot:
    α_true = 1 + |s|    (NOT |s| itself)

MLE works directly on raw data, not on a plot. It is unbiased and gives
correct α regardless of how you later choose to visualise the result.

Usage
-----
    python bpl_mle.py                         # uses default path below
    python bpl_mle.py my_data.csv             # pass your own CSV
    python bpl_mle.py my_data.csv intensity   # specify column name
    python bpl_mle.py my_data.csv intensity 0.5   # also specify x_min
"""

import sys, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import minimize, differential_evolution, minimize_scalar
from scipy.stats import kstest

warnings.filterwarnings("ignore")

# ── User settings (override via command-line args) ────────────────────────────
DATA_PATH  = "/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Comapre_flares/pre_flare_intensity_distribution/Flux_distrib/full_data_unsorted.csv"
COL_NAME   = "intensity"
X_MIN      = 0.33          # events below this are excluded (not real flares)
OUT_PNG    = "broken_powerlaw_mle.png"

if len(sys.argv) > 1: DATA_PATH = sys.argv[1]
if len(sys.argv) > 2: COL_NAME  = sys.argv[2]
if len(sys.argv) > 3: X_MIN     = float(sys.argv[3])


# ══════════════════════════════════════════════════════════════════════════════
# 1. LOAD DATA
# ══════════════════════════════════════════════════════════════════════════════
df       = pd.read_csv(DATA_PATH)
data_all = df[COL_NAME].dropna().values
data     = np.sort(data_all[data_all >= X_MIN])
n        = len(data)
x_max    = data.max()

print(f"\n{'='*55}")
print(f"  Dataset         : {DATA_PATH}")
print(f"  Column          : {COL_NAME}")
print(f"  Total events    : {len(data_all)}")
print(f"  Events ≥ {X_MIN}  : {n}")
print(f"  Range (used)    : [{data.min():.4f}, {x_max:.4f}]")
print(f"{'='*55}\n")


# ══════════════════════════════════════════════════════════════════════════════
# 2. BROKEN POWER LAW MODEL
# ══════════════════════════════════════════════════════════════════════════════
#
# PDF:   f(x) = C · x^{-α1}                     x_min ≤ x < x_b
#        f(x) = C · x_b^{α2-α1} · x^{-α2}       x ≥ x_b
#
# Normalisation Z = ∫ f(x)/C dx  over [x_min, ∞)
#
#   Z = ∫_{x_min}^{x_b} x^{-α1} dx   +   ∫_{x_b}^{∞} x_b^{α2-α1} x^{-α2} dx
#
#   Piece 1:  [x^{1-α1} / (1-α1)]_{x_min}^{x_b}
#             = (x_b^{1-α1} - x_min^{1-α1}) / (1-α1)      [α1 ≠ 1]
#             = ln(x_b / x_min)                              [α1 = 1 exactly]
#
#   Piece 2:  x_b^{α2-α1} · [x^{1-α2}/(1-α2)]_{x_b}^{∞}
#             = x_b^{1-α1} / (α2-1)                         [requires α2 > 1]
#
# IMPORTANT: the denominator in Piece 1 is (1-α1), NOT (α1-1).
# Swapping the sign is a common bug that silently inflates the likelihood.

def normalisation_Z(a1, a2, xb):
    """Normalisation integral Z.  α1 unrestricted; requires α2 > 1."""
    if abs(a1 - 1.0) < 1e-9:
        piece1 = np.log(xb / X_MIN)
    else:
        piece1 = (xb**(1-a1) - X_MIN**(1-a1)) / (1 - a1)   # ← correct sign
    piece2 = xb**(1-a1) / (a2 - 1)
    return piece1 + piece2


def neg_log_lik(params):
    """Negative log-likelihood for minimisation."""
    a1, a2, xb = params
    # α2 > 1 ensures upper integral converges; xb must be inside data range
    if a2 <= 1.001 or xb <= X_MIN + 0.005 or xb >= x_max - 0.005:
        return 1e15
    Z = normalisation_Z(a1, a2, xb)
    if Z <= 0 or not np.isfinite(Z):
        return 1e15

    lo = data[data < xb]
    hi = data[data >= xb]

    # log p(xi): cancel C by absorbing -log(Z) as n*(-log Z)
    ll  = -n * np.log(Z)
    if len(lo): ll += np.sum(-a1 * np.log(lo))
    if len(hi): ll += np.sum((a2 - a1)*np.log(xb) - a2*np.log(hi))
    return -ll


# ══════════════════════════════════════════════════════════════════════════════
# 3. OPTIMISATION  (global DE search → local Nelder-Mead polish)
# ══════════════════════════════════════════════════════════════════════════════
print("Fitting broken power law via MLE ...")

bounds = [
    (1.0, 2.0),                          # α1: no lower constraint (data may rise)
    (2.001, 5.0),                        # α2: must be > 1
    (0.8, 1.1),        # x_b
]

res_de = differential_evolution(
    neg_log_lik, bounds,
    seed=42, maxiter=6000, tol=1e-13,
    popsize=30, mutation=(0.5, 1.5), recombination=0.7,
    workers=1,
)

res = minimize(
    neg_log_lik, res_de.x, method="Nelder-Mead",
    options={"xatol": 1e-12, "fatol": 1e-12, "maxiter": 500_000},
)

a1_mle, a2_mle, xb_mle = res.x
max_ll = -res.fun
Z_mle  = normalisation_Z(a1_mle, a2_mle, xb_mle)


# ══════════════════════════════════════════════════════════════════════════════
# 4. UNCERTAINTY ESTIMATES — numerical Hessian → Fisher information matrix
# ══════════════════════════════════════════════════════════════════════════════
def num_hessian(f, x0, eps=5e-4):
    d = len(x0); H = np.zeros((d, d))
    for i in range(d):
        for j in range(d):
            pp=x0.copy(); pp[i]+=eps; pp[j]+=eps
            pm=x0.copy(); pm[i]+=eps; pm[j]-=eps
            mp=x0.copy(); mp[i]-=eps; mp[j]+=eps
            mm=x0.copy(); mm[i]-=eps; mm[j]-=eps
            H[i,j] = (f(pp)-f(pm)-f(mp)+f(mm)) / (4*eps**2)
    return H

try:
    H    = num_hessian(neg_log_lik, res.x)
    cov  = np.linalg.inv(H)
    se   = np.sqrt(np.abs(np.diag(cov)))
    ci95 = 1.96 * se
except Exception as exc:
    print(f"  Warning: Hessian failed ({exc}); SE set to NaN")
    se = ci95 = np.array([np.nan, np.nan, np.nan])


# ══════════════════════════════════════════════════════════════════════════════
# 5. SIMPLE POWER LAW COMPARISON (AIC / BIC)
# ══════════════════════════════════════════════════════════════════════════════
def neg_ll_simple(a):
    if a <= 1.001: return 1e15
    return -(n*np.log(a-1) + n*(a-1)*np.log(X_MIN) - a*np.sum(np.log(data)))

res_s  = minimize_scalar(neg_ll_simple, bounds=(1.01, 12), method='bounded')
a_spl  = res_s.x
ll_spl = -res_s.fun

aic_bpl = 2*3 - 2*max_ll;    bic_bpl = 3*np.log(n) - 2*max_ll
aic_spl = 2*1 - 2*ll_spl;    bic_spl = 1*np.log(n) - 2*ll_spl


# ══════════════════════════════════════════════════════════════════════════════
# 6. PROFILE LIKELIHOOD FOR x_b (more reliable CI than Hessian for this param)
# ══════════════════════════════════════════════════════════════════════════════
xb_grid   = np.linspace(X_MIN + 0.08, x_max - 0.08, 100)
profile   = np.full(len(xb_grid), np.nan)

for k, xb_try in enumerate(xb_grid):
    def _nll(ab):
        return neg_log_lik([ab[0], ab[1], xb_try])
    r = minimize(_nll, [a1_mle, a2_mle], method='Nelder-Mead',
                 options={'xatol':1e-9,'fatol':1e-9,'maxiter':50000,'disp':False})
    if r.fun < 1e14:
        profile[k] = -r.fun

# 95% CI from chi^2 cutoff: profile >= max_ll - 1.92  (chi^2_1 / 2)
valid_ci = xb_grid[profile >= max_ll - 1.92]
xb_ci_lo = valid_ci.min() if len(valid_ci) else np.nan
xb_ci_hi = valid_ci.max() if len(valid_ci) else np.nan


# ══════════════════════════════════════════════════════════════════════════════
# 7. GOODNESS OF FIT  — KS test on raw data
# ══════════════════════════════════════════════════════════════════════════════
def bpl_cdf(x, a1, a2, xb):
    Z   = normalisation_Z(a1, a2, xb)
    cdf = np.zeros_like(x, dtype=float)
    lo  = (x >= X_MIN) & (x < xb)
    hi  = x >= xb
    if lo.any():
        if abs(a1-1) < 1e-9:
            cdf[lo] = np.log(x[lo]/X_MIN) / Z
        else:
            cdf[lo] = (x[lo]**(1-a1) - X_MIN**(1-a1)) / ((1-a1) * Z)
    if hi.any():
        I_lo = (np.log(xb/X_MIN) if abs(a1-1)<1e-9
                else (xb**(1-a1) - X_MIN**(1-a1)) / (1-a1))
        # ∫_{xb}^{x} x_b^{a2-a1} t^{-a2} dt = x_b^{a2-a1}·[t^{1-a2}/(1-a2)]_{xb}^{x}
        I_hi = xb**(a2-a1) * (xb**(1-a2) - x[hi]**(1-a2)) / (1-a2)
        cdf[hi] = (I_lo + I_hi) / Z
    return np.clip(cdf, 0.0, 1.0)

ks_D, ks_p = kstest(data, lambda x: bpl_cdf(x, a1_mle, a2_mle, xb_mle))


# ══════════════════════════════════════════════════════════════════════════════
# 8. PRINT RESULTS
# ══════════════════════════════════════════════════════════════════════════════
#
#  READING SLOPES OFF PLOTS:
#  ─────────────────────────
#  PDF log-log slope  = -α      → α = |slope|
#  CCDF log-log slope = -(α-1)  → α = 1 + |slope|
#
#  So below x_b, the CCDF slope you see ≈ -(α1-1).
#  Above x_b, the CCDF slope you see ≈ -(α2-1).
#  The PDF slopes are one unit steeper.
#
#  For our fit:  α1=0.821  → CCDF slope below break ≈ +(0.179)  [slight uptick]
#                α2=5.876  → CCDF slope above break ≈ -4.876

print(f"{'='*55}")
print(f"  RESULTS")
print(f"{'='*55}")
print(f"\n  Broken Power Law (MLE):")
print(f"    α₁  (PDF slope below break) = {a1_mle:.3f} ± {se[0]:.3f}")
print(f"    α₂  (PDF slope above break) = {a2_mle:.3f} ± {se[1]:.3f}")
print(f"    x_b (break point)           = {xb_mle:.3f} ± {se[2]:.3f}")
print(f"    x_b 95% CI (profile)        = [{xb_ci_lo:.3f}, {xb_ci_hi:.3f}]")
print(f"    Log-likelihood              = {max_ll:.3f}")
print(f"    AIC / BIC                   = {aic_bpl:.2f} / {bic_bpl:.2f}")
print(f"\n  What you'd READ off a CCDF plot:")
print(f"    Slope below x_b  ≈ {-(a1_mle-1):+.3f}  (α₁-1 = {a1_mle-1:.3f})")
print(f"    Slope above x_b  ≈ {-(a2_mle-1):+.3f}  (α₂-1 = {a2_mle-1:.3f})")
print(f"\n  Simple Power Law (MLE):")
print(f"    α                           = {a_spl:.3f}")
print(f"    AIC / BIC                   = {aic_spl:.2f} / {bic_spl:.2f}")
print(f"\n  Model selection:")
print(f"    ΔAIC = {aic_spl-aic_bpl:.1f}  |  ΔBIC = {bic_spl-bic_bpl:.1f}")
pref = "Broken PL strongly preferred" if aic_spl - aic_bpl > 10 else \
       "Broken PL preferred" if aic_spl - aic_bpl > 2 else "No strong preference"
print(f"    → {pref}")
print(f"\n  KS test: D = {ks_D:.4f},  p = {ks_p:.4f}")
print(f"{'='*55}\n")


# ══════════════════════════════════════════════════════════════════════════════
# 9. FIGURE
# ══════════════════════════════════════════════════════════════════════════════
BG="#0d0f18"; PL="#131625"; CG="#252a40"; CT="#dce3f0"
C_dat="#7ec8e3"; C_bpl="#f4a261"; C_spl="#e63946"; C_brk="#a8dadc"

fig = plt.figure(figsize=(16, 10))
fig.patch.set_facecolor(BG)
gs  = gridspec.GridSpec(2, 3, figure=fig,
                        hspace=0.46, wspace=0.36,
                        left=0.07, right=0.97, top=0.91, bottom=0.09)

def _ax(ax):
    ax.set_facecolor(PL)
    ax.tick_params(colors=CT, labelsize=8.5)
    ax.xaxis.label.set_color(CT); ax.yaxis.label.set_color(CT)
    ax.title.set_color(CT)
    for sp in ax.spines.values(): sp.set_color(CG)
    ax.grid(True, color=CG, lw=0.5, ls="--", alpha=0.6)
    return ax

x_lo  = X_MIN * 0.90
x_hi  = x_max * 1.08
x_plt = np.logspace(np.log10(X_MIN), np.log10(x_max * 1.04), 800)

C_hat = 1.0 / Z_mle
def bpl_pdf_fn(x):
    pdf = np.zeros_like(x, dtype=float)
    m   = x < xb_mle
    pdf[ m] = C_hat * x[m]**(-a1_mle)
    pdf[~m] = C_hat * xb_mle**(a2_mle-a1_mle) * x[~m]**(-a2_mle)
    return pdf

C_s = (a_spl - 1) * X_MIN**(a_spl - 1)
def spl_pdf_fn(x): return C_s * x**(-a_spl)

# ── A: Log-log PDF ────────────────────────────────────────────────────────────
ax1 = _ax(fig.add_subplot(gs[0, :2]))

log_bins  = np.logspace(np.log10(X_MIN), np.log10(x_max), 18)
cnt, ed   = np.histogram(data, bins=log_bins)
bw        = ed[1:] - ed[:-1]
bc        = np.sqrt(ed[:-1]*ed[1:])
pdf_h     = cnt / (n * bw)
v         = cnt > 0

ax1.errorbar(bc[v], pdf_h[v], yerr=np.sqrt(cnt[v])/(n*bw[v]),
             fmt='o', color=C_dat, ms=6, elinewidth=1.3, capsize=4,
             label=f"Data (log-binned PDF, n={n})", zorder=5)

ax1.loglog(x_plt, bpl_pdf_fn(x_plt), color=C_bpl, lw=2.5, zorder=6,
           label=(f"Broken PL:  α₁={a1_mle:.2f}±{se[0]:.2f},  "
                  f"α₂={a2_mle:.2f}±{se[1]:.2f},  x_b={xb_mle:.2f}"))
ax1.loglog(x_plt, spl_pdf_fn(x_plt), color=C_spl, lw=1.6, ls="--", zorder=4,
           label=f"Simple PL:  α={a_spl:.2f}")

ax1.axvline(xb_mle, color=C_brk, lw=2, ls=":", alpha=0.9,
            label=f"Break  x_b={xb_mle:.3f}")
ax1.axvline(X_MIN, color="white", lw=1.0, ls=":", alpha=0.3)

# Annotate PDF slopes
yrange = ax1.get_ylim()
mid_lo = np.sqrt(X_MIN * xb_mle)
mid_hi = np.sqrt(xb_mle * x_max)
ax1.text(mid_lo, bpl_pdf_fn(np.array([mid_lo]))[0]*2.8,
         f"PDF slope\n= −{a1_mle:.2f}", color=C_bpl, fontsize=7.5,
         ha='center', va='bottom')
ax1.text(mid_hi, bpl_pdf_fn(np.array([mid_hi]))[0]*2.8,
         f"PDF slope\n= −{a2_mle:.2f}", color=C_bpl, fontsize=7.5,
         ha='center', va='bottom')

ax1.set_xlabel("Peak Flux Intensity  (erg)", fontsize=9.5)
ax1.set_ylabel("PDF  f(x)", fontsize=9.5)
ax1.set_title("Log-log PDF  (slopes here = −α)", fontsize=11, pad=8)
ax1.legend(fontsize=8.5, framealpha=0.3, facecolor=BG, labelcolor=CT, loc="upper right")

# ── B: Log-log CCDF ──────────────────────────────────────────────────────────
ax2 = _ax(fig.add_subplot(gs[0, 2]))

x_emp   = np.sort(data)
ccdf_emp = 1.0 - np.arange(1, n+1)/n

ccdf_bpl = 1.0 - bpl_cdf(x_emp, a1_mle, a2_mle, xb_mle)
ccdf_spl = (x_emp / X_MIN)**(-(a_spl-1))

ax2.loglog(x_emp, np.maximum(ccdf_emp, 1e-4), color=C_dat, lw=0,
           marker='.', ms=4.5, alpha=0.8, label="Empirical CCDF")
ax2.loglog(x_emp, np.maximum(ccdf_bpl, 1e-4), color=C_bpl, lw=2.3,
           label="Broken PL")
ax2.loglog(x_emp, np.maximum(ccdf_spl, 1e-4), color=C_spl, lw=1.5, ls="--",
           label="Simple PL")
ax2.axvline(xb_mle, color=C_brk, lw=1.5, ls=":", alpha=0.8)

# Annotate CCDF slopes (= -(α-1))
x_ann_lo = np.sqrt(X_MIN * xb_mle)
x_ann_hi = np.sqrt(xb_mle * x_max)
idx_lo = np.searchsorted(x_emp, x_ann_lo)
idx_hi = np.searchsorted(x_emp, x_ann_hi)

if 0 < idx_lo < len(ccdf_bpl) and ccdf_bpl[idx_lo] > 0:
    ax2.text(x_ann_lo, ccdf_bpl[idx_lo]*3,
             f"CCDF slope\n≈{-(a1_mle-1):+.2f}", color=C_bpl,
             fontsize=7, ha='center', va='bottom')
if 0 < idx_hi < len(ccdf_bpl) and ccdf_bpl[idx_hi] > 0:
    ax2.text(x_ann_hi, ccdf_bpl[idx_hi]*3,
             f"CCDF slope\n≈{-(a2_mle-1):+.2f}", color=C_bpl,
             fontsize=7, ha='center', va='bottom')

ax2.set_xlabel("Peak Flux  (erg)", fontsize=9)
ax2.set_ylabel("P(X > x)  [CCDF]", fontsize=9)
ax2.set_title("Log-log CCDF  (slopes here = −(α−1))", fontsize=9.5, pad=8)
ax2.legend(fontsize=8, framealpha=0.3, facecolor=BG, labelcolor=CT)

# ── C: Raw histogram with threshold ──────────────────────────────────────────
ax3 = _ax(fig.add_subplot(gs[1, 0]))

bins_lin   = np.linspace(0, x_max*1.05, 30)
cnt_all, ed_all = np.histogram(data_all, bins=bins_lin)
bw_all = ed_all[1:] - ed_all[:-1]
bc_all = 0.5*(ed_all[:-1]+ed_all[1:])

ax3.bar(bc_all, cnt_all, width=bw_all*0.85, color=C_dat, alpha=0.5,
        label="All events")
ax3.bar(bc_all[bc_all >= X_MIN], cnt_all[bc_all >= X_MIN],
        width=bw_all[bc_all >= X_MIN]*0.85, color=C_bpl, alpha=0.7,
        label=f"Used (≥ {X_MIN})")
ax3.axvline(X_MIN, color="white", lw=1.5, ls="--", alpha=0.7,
            label=f"x_min={X_MIN}")
ax3.axvline(xb_mle, color=C_brk, lw=1.8, ls=":", alpha=0.9,
            label=f"x_b={xb_mle:.2f}")
ax3.set_xlabel("Peak Flux  (erg)", fontsize=9)
ax3.set_ylabel("Count", fontsize=9)
ax3.set_title("Linear histogram (all data)", fontsize=9.5)
ax3.legend(fontsize=8, framealpha=0.3, facecolor=BG, labelcolor=CT)

# ── D: Profile likelihood for x_b ────────────────────────────────────────────
ax4 = _ax(fig.add_subplot(gs[1, 1]))

valid_mask = np.isfinite(profile)
ax4.plot(xb_grid[valid_mask], profile[valid_mask], color=C_bpl, lw=2)
ax4.axvline(xb_mle,    color=C_brk, lw=2, ls=":",  alpha=0.9,
            label=f"MLE  x_b={xb_mle:.3f}")
ax4.axhline(max_ll - 1.92, color="white", lw=1.0, ls="--", alpha=0.5,
            label="95% CI boundary")
ax4.fill_between(xb_grid[valid_mask], profile[valid_mask],
                 max_ll - 1.92,
                 where=(profile[valid_mask] >= max_ll - 1.92),
                 color=C_bpl, alpha=0.18)
ax4.axvline(xb_ci_lo, color="white", lw=0.8, ls=":", alpha=0.4)
ax4.axvline(xb_ci_hi, color="white", lw=0.8, ls=":", alpha=0.4)
ax4.set_xlabel("Break point  x_b", fontsize=9)
ax4.set_ylabel("Profile log-likelihood", fontsize=9)
ax4.set_title(f"Profile LL for x_b  [95% CI: {xb_ci_lo:.2f}–{xb_ci_hi:.2f}]",
              fontsize=9.5)
ax4.legend(fontsize=8, framealpha=0.3, facecolor=BG, labelcolor=CT)

# ── E: Results + explanation box ─────────────────────────────────────────────
ax5 = fig.add_subplot(gs[1, 2])
ax5.set_facecolor(PL); ax5.set_xticks([]); ax5.set_yticks([])
for sp in ax5.spines.values(): sp.set_color(CG)
ax5.set_title("Results & Slope Guide", fontsize=9.5, color=CT, pad=5)

txt = [
    f"n={n}   x_min={X_MIN}",
    "",
    "── Broken PL (MLE) ──────────",
    f"  α₁  = {a1_mle:.3f} ± {se[0]:.3f}",
    f"  α₂  = {a2_mle:.3f} ± {se[1]:.3f}",
    f"  x_b = {xb_mle:.3f}  [{xb_ci_lo:.2f}–{xb_ci_hi:.2f}]",
    f"  AIC={aic_bpl:.1f}  BIC={bic_bpl:.1f}",
    "",
    "── Simple PL (MLE) ──────────",
    f"  α   = {a_spl:.3f}",
    f"  AIC={aic_spl:.1f}  BIC={bic_spl:.1f}",
    "",
    f"ΔAIC={aic_spl-aic_bpl:.1f}  ΔBIC={bic_spl-bic_bpl:.1f}",
    f"KS: D={ks_D:.3f}  p={ks_p:.3f}",
    "",
    "── Slope cheatsheet ─────────",
    "  PDF  slope = -α",
    "  CCDF slope = -(α-1)",
    "  → α = 1 + |CCDF slope|",
    "",
    f"  Below break (CCDF slope):",
    f"    ≈ {-(a1_mle-1):+.3f}  (α₁-1={a1_mle-1:.3f})",
    f"  Above break (CCDF slope):",
    f"    ≈ {-(a2_mle-1):+.3f}  (α₂-1={a2_mle-1:.3f})",
]
ax5.text(0.04, 0.97, "\n".join(txt), transform=ax5.transAxes,
         fontsize=7.6, va="top", family="monospace", color=CT, linespacing=1.55)

plt.suptitle(
    "Pre-Flare Transient Peak Flux — Broken Power Law (MLE)\n"
    "PDF slope = −α  |  CCDF slope = −(α−1)",
    fontsize=12, color=CT, y=0.97)

fig.savefig(OUT_PNG, dpi=160, bbox_inches="tight", facecolor=BG)
plt.close()
print(f"Figure saved → {OUT_PNG}")
