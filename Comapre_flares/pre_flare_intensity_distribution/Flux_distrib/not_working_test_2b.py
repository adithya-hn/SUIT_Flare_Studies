import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution, curve_fit
from scipy.stats import kstest

# ══════════════════════════════════════════════════════
# 0. LOAD DATA  — replace with your actual array
# ══════════════════════════════════════════════════════
# data = np.loadtxt("your_energies.txt")
np.random.seed(42)
data = np.concatenate([
    np.random.pareto(0.2, 30) * 0.015 + 0.01,   # flat
    np.random.pareto(1.5, 70) * 0.04  + 0.08,   # mid
    np.random.pareto(5.0, 20) * 0.08  + 0.70,   # steep
])
data = np.sort(data[data > 0])

# ══════════════════════════════════════════════════════
# 1. EMPIRICAL CCDF
# ══════════════════════════════════════════════════════
n    = len(data)
ccdf = 1.0 - np.arange(1, n + 1) / n
lx   = np.log10(data)
ly   = np.log10(np.clip(ccdf, 1e-10, None))

print(f"N = {n} events")
print(f"Energy range: [{data.min():.4f}, {data.max():.4f}]")

# ══════════════════════════════════════════════════════
# 2. PIECEWISE MODELS  (continuous in log-log space)
# ══════════════════════════════════════════════════════
def seg2_model(lx, lbp, a1, a2, lC):
    """1 break → 2 power laws"""
    y = np.empty_like(lx, dtype=float)
    m1, m2 = lx <= lbp, lx > lbp
    y[m1]  = lC + a1 * (lx[m1] - lx[0])
    y_bp   = lC + a1 * (lbp    - lx[0])
    y[m2]  = y_bp + a2 * (lx[m2] - lbp)
    return y

def seg3_model(lx, lbp1, lbp2, a1, a2, a3, lC):
    """2 breaks → 3 power laws"""
    y = np.empty_like(lx, dtype=float)
    m1 = lx <= lbp1
    m2 = (lx > lbp1) & (lx <= lbp2)
    m3 = lx > lbp2
    y[m1]  = lC + a1 * (lx[m1] - lx[0])
    y_bp1  = lC + a1 * (lbp1 - lx[0])
    y[m2]  = y_bp1 + a2 * (lx[m2] - lbp1)
    y_bp2  = y_bp1 + a2 * (lbp2 - lbp1)
    y[m3]  = y_bp2 + a3 * (lx[m3] - lbp2)
    return y

# ══════════════════════════════════════════════════════
# 3. GLOBAL OPTIMIZER  (differential evolution)
# ══════════════════════════════════════════════════════
lx0, lx1 = lx.min(), lx.max()
span      = lx1 - lx0

def rss(y_pred):
    return np.sum((ly - y_pred)**2)

# --- 1-segment baseline ---
coeffs   = np.polyfit(lx, ly, 1)
y_1seg   = np.polyval(coeffs, lx)
rss_1    = rss(y_1seg)
alpha_1  = -coeffs[0]

# --- 2-segment fit ---
def cost_2seg(p):
    lbp, a1, a2, lC = p
    return rss(seg2_model(lx, lbp, a1, a2, lC))

res2 = differential_evolution(
    cost_2seg,
    bounds=[(lx0+0.1*span, lx1-0.1*span),
            (-3, 2), (-10, 0), (ly.min(), ly.max())],
    seed=42, maxiter=3000, tol=1e-12, popsize=25, polish=True
)
y_2seg = seg2_model(lx, *res2.x)
rss_2  = rss(y_2seg)

# --- 3-segment fit ---
def cost_3seg(p):
    lbp1, lbp2, a1, a2, a3, lC = p
    if lbp1 >= lbp2:
        return 1e10
    return rss(seg3_model(lx, lbp1, lbp2, a1, a2, a3, lC))

res3 = differential_evolution(
    cost_3seg,
    bounds=[(lx0+0.05*span, lx0+0.45*span),   # lbp1
            (lx0+0.40*span, lx1-0.05*span),    # lbp2
            (-2,  2),   # a1 flat
            (-5,  0),   # a2 moderate
            (-15, 0),   # a3 steep
            (ly.min(), ly.max())],
    seed=42, maxiter=5000, tol=1e-12, popsize=30, polish=True
)
y_3seg = seg3_model(lx, *res3.x)
rss_3  = rss(y_3seg)

# ══════════════════════════════════════════════════════
# 4. MODEL SELECTION  — BIC & AIC
#    Lower BIC/AIC = better model
# ══════════════════════════════════════════════════════
def bic_aic(rss_val, k, n):
    """k = number of free parameters"""
    sigma2 = rss_val / n
    ll     = -0.5 * n * np.log(2 * np.pi * sigma2) - rss_val / (2 * sigma2)
    bic    = -2 * ll + k * np.log(n)
    aic    = -2 * ll + 2 * k
    aicc   = aic + (2*k*(k+1)) / (n - k - 1)   # corrected AIC for small n
    return bic, aic, aicc

bic1, aic1, aicc1 = bic_aic(rss_1, 2, n)
bic2, aic2, aicc2 = bic_aic(rss_2, 4, n)
bic3, aic3, aicc3 = bic_aic(rss_3, 6, n)

print("\n╔══════════════════════════════════════════════════╗")
print("║          MODEL SELECTION SUMMARY                ║")
print("╠══════════════════════════════════════════════════╣")
print(f"║  Model     k    BIC       AIC      AICc        ║")
print(f"║  1-segment 2  {bic1:8.2f}  {aic1:8.2f}  {aicc1:8.2f}      ║")
print(f"║  2-segment 4  {bic2:8.2f}  {aic2:8.2f}  {aicc2:8.2f}      ║")
print(f"║  3-segment 6  {bic3:8.2f}  {aic3:8.2f}  {aicc3:8.2f}      ║")
print("╚══════════════════════════════════════════════════╝")

bics   = [bic1, bic2, bic3]
best_k = np.argmin(bics) + 1
delta  = [b - min(bics) for b in bics]
print(f"\n  ✔ Best model: {best_k} segment(s)")
print(f"  ΔBIC vs best:  1-seg={delta[0]:.1f},  "
      f"2-seg={delta[1]:.1f},  3-seg={delta[2]:.1f}")
print("  (ΔBIC > 10 = strong evidence against that model)")

# ══════════════════════════════════════════════════════
# 5. PRINT BEST-FIT PARAMETERS
# ══════════════════════════════════════════════════════
lbp1_B, lbp2_B, a1_B, a2_B, a3_B, lC_B = res3.x
lbp_2,  a1_2,   a2_2, lC_2             = res2.x

print("\n─── 2-segment parameters ───")
print(f"  Break energy   : {10**lbp_2:.4f}  erg cm⁻²")
print(f"  α₁ (below)     : {-a1_2:.3f}")
print(f"  α₂ (above)     : {-a2_2:.3f}")

print("\n─── 3-segment parameters ───")
print(f"  Break energy 1 : {10**lbp1_B:.4f}  erg cm⁻²")
print(f"  Break energy 2 : {10**lbp2_B:.4f}  erg cm⁻²")
print(f"  α₁ (flat)      : {-a1_B:.3f}")
print(f"  α₂ (mid)       : {-a2_B:.3f}")
print(f"  α₃ (steep)     : {-a3_B:.3f}")

# ══════════════════════════════════════════════════════
# 6. KS TEST — goodness of fit for best model
# ══════════════════════════════════════════════════════
y_best = y_3seg if best_k == 3 else y_2seg
resid  = ly - y_best
ks_stat, ks_p = kstest(resid, 'norm',
                        args=(resid.mean(), resid.std()))
print(f"\n  KS test on residuals: stat={ks_stat:.3f}, p={ks_p:.3f}")
print(f"  {'✔ Residuals look normal' if ks_p > 0.05 else '⚠ Residuals non-normal — check fit'}")

# ══════════════════════════════════════════════════════
# 7. PLOT — all models + residuals
# ══════════════════════════════════════════════════════
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 9),
                                gridspec_kw={'height_ratios': [3, 1]},
                                sharex=False)

lx_fine = np.linspace(lx.min(), lx.max(), 2000)

# ── Main panel ──
ax1.plot(data, ccdf, 'k+', ms=5, zorder=5, label='SUIT transient events')

# 2-segment
y2f = seg2_model(lx_fine, *res2.x)
ax1.plot(10**lx_fine, 10**y2f, 'b--', lw=1.8,
         label=fr'2-seg  $\alpha_1={-a1_2:.2f},\ \alpha_2={-a2_2:.2f}$')

# 3-segment
y3f = seg3_model(lx_fine, *res3.x)
ax1.plot(10**lx_fine, 10**y3f, 'r-', lw=2,
         label=fr'3-seg  $\alpha_1={-a1_B:.2f},\ \alpha_2={-a2_B:.2f},\ \alpha_3={-a3_B:.2f}$')

# Breakpoints
ax1.axvline(10**lbp1_B, color='tomato',   ls=':', lw=1.5,
            label=fr'$E_{{b1}}={10**lbp1_B:.3f}$')
ax1.axvline(10**lbp2_B, color='darkorange', ls=':', lw=1.5,
            label=fr'$E_{{b2}}={10**lbp2_B:.3f}$')

# BIC annotation
best_label = {1:"1-seg", 2:"2-seg", 3:"3-seg"}[best_k]
ax1.text(0.03, 0.05,
         f"Best model (BIC): {best_label}\n"
         f"ΔBIC(2-seg)={delta[1]:.1f}  ΔBIC(3-seg)={delta[2]:.1f}",
         transform=ax1.transAxes, fontsize=10,
         bbox=dict(boxstyle='round', fc='wheat', alpha=0.7))

ax1.set_xscale('log'); ax1.set_yscale('log')
ax1.set_ylabel('CCDF', fontsize=13)
ax1.set_title('CCDF of pre-flare transients — automatic model selection', fontsize=13)
ax1.legend(fontsize=9, loc='lower left')
ax1.grid(True, which='both', ls='--', alpha=0.35)

# ── Residual panel ──
ax2.axhline(0, color='k', lw=0.8)
ax2.plot(data, ly - y_3seg, 'r.', ms=4, alpha=0.6, label='3-seg residuals')
ax2.plot(data, ly - y_2seg, 'b.', ms=4, alpha=0.4, label='2-seg residuals')
ax2.set_xscale('log')
ax2.set_xlabel(r'Energy ($erg^{-1} \times cm^{-2}$)', fontsize=13)
ax2.set_ylabel('Residual\n(log space)', fontsize=10)
ax2.legend(fontsize=9)
ax2.grid(True, which='both', ls='--', alpha=0.35)

plt.tight_layout()
plt.savefig('ccdf_model_selection.png', dpi=150)
plt.show()