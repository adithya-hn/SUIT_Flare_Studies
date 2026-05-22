import numpy as np
import matplotlib.pyplot as plt
import scienceplots
import matplotlib.dates as mdates
plt.style.use('science')

# ---------------------------------------------------
# YOUR DATA

# ---------------------------------------------------
data = np.loadtxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Comapre_flares/pre_flare_intensity_distribution/Flux_distrib/Peak_Flux.csv')

import numpy as np
from scipy.optimize import minimize_scalar
data_=np.sort(data)
data = np.sort(data[data > 0.36])

EA   = 0.36          # fixed theoretical lower bound — NOT data.min()
xmax = data.max()

def logL_truncated(alpha, x, xlo, xhi):
    if alpha <= 1:
        return -np.inf
    N    = len(x)
    norm = (xhi**(1 - alpha) - xlo**(1 - alpha)) / (1 - alpha)
    if norm <= 0 or not np.isfinite(norm):
        return -np.inf
    return -alpha * np.sum(np.log(x)) - N * np.log(norm)

def mle_alpha_truncated(x, xlo, xhi):
    """Exact MLE via scipy — no grid quantisation."""
    result = minimize_scalar(
        lambda a: -logL_truncated(a, x, xlo, xhi),
        bounds=(1.001, 8.0),
        method='bounded'
    )
    return result.x

# Break point search over interior data points
candidate_breaks = data[10:-10]

best_ll     = -np.inf
best_params = None

for xb in candidate_breaks:
    seg1 = data[(data >= EA)   & (data < xb)]
    seg2 = data[(data >= xb)   & (data <= xmax)]

    if len(seg1) < 10 or len(seg2) < 10:
        continue

    a1 = mle_alpha_truncated(seg1, EA, xb)
    a2 = mle_alpha_truncated(seg2, xb, xmax)

    ll = (logL_truncated(a1, seg1, EA, xb)
        + logL_truncated(a2, seg2, xb, xmax))

    if ll > best_ll:
        best_ll     = ll
        best_params = (xb, a1, a2)

xbest, alpha1, alpha2 = best_params

print(f"Best break  = {xbest:.5f}")
print(f"Alpha1      = {alpha1:.5f}  (n1 = {np.sum(data < xbest)})")
print(f"Alpha2      = {alpha2:.5f}  (n2 = {np.sum(data >= xbest)})")
print(f"LogL        = {best_ll:.5f}")
print()

def empirical_ccdf(E):
    n = len(E)
    return np.arange(n, 0, -1) / n


import numpy as np
import matplotlib.pyplot as plt

# ── Inputs (from your MLE fit) ───────────────────────────
EA, EB   = 0.36, data.max()
N        = len(data)
xb       = 1.03611          # single break point
a1, a2   = 1.695, 2.983     # MLE indices
a1_lo, a1_hi = 1.131, 2.240 # 68% bootstrap CI
a2_lo, a2_hi = 2.391, 3.690

seg1 = data[data < xb];  n1 = len(seg1)
seg2 = data[data >= xb]; n2 = len(seg2)

# ── Correct truncated PL CCDF (scaled to N fractions) ────
def ccdf_trunc(x, a, xlo, xhi, scale):
    norm = (xhi**(1-a) - xlo**(1-a)) / (1-a)
    return scale * np.clip(
        1 - (x**(1-a) - xlo**(1-a)) / ((1-a)*norm), 0, 1)

x1 = np.logspace(np.log10(EA),  np.log10(xb), 400)
x2 = np.logspace(np.log10(xb),  np.log10(EB), 400)
c1 = ccdf_trunc(x1, a1, EA, xb, n1/N)
c2 = ccdf_trunc(x2, a2, xb, EB, n2/N)

# # 68% bands (recompute from bootstrap arrays b1, b2)
# lo1 = np.percentile([ccdf_trunc(x1,b,EA,xb,n1/N) for b in b1], 16, axis=0)
# hi1 = np.percentile([ccdf_trunc(x1,b,EA,xb,n1/N) for b in b1], 84, axis=0)
# lo2 = np.percentile([ccdf_trunc(x2,b,xb,EB,n2/N) for b in b2], 16, axis=0)
# hi2 = np.percentile([ccdf_trunc(x2,b,xb,EB,n2/N) for b in b2], 84, axis=0)

# ── Plot ──────────────────────────────────────────────────


# fig, ax = plt.subplots(figsize=(7, 5))

# x_emp    = np.sort(data)
# ccdf_emp = np.arange(N, 0, -1) / N
# ax.plot(x_emp, ccdf_emp, '+', color='#2C2C2A', ms=5, mew=0.9,zorder=5, label='SUIT Mg\u202fII h transients')

# # ax.fill_between(x1, lo1, hi1, color='#185FA5', alpha=0.15, linewidth=0)
# # ax.fill_between(x2, lo2, hi2, color='#D85A30', alpha=0.15, linewidth=0)

# ax.plot(x1, c1, color='#185FA5', lw=2.0, zorder=4, label=rf'$\alpha_1={a1:.2f}$'  r'  ($x < x_b$)')
# ax.plot(x2, c2, color='#D85A30', lw=2.0, zorder=4, label=rf'$\alpha_2={a2:.2f}$'r'  ($x \geq x_b$)')

# ax.axvline(xb, color='#BA7517', lw=1.4, ls='--', zorder=3,
#            label=rf'$x_b = {xb:.3f}$ erg s$^{{-1}}$ cm$^{{-2}}$')

# ax.set_xscale('log'); ax.set_yscale('log')
# ax.set_xlim(0.28, 4.0); ax.set_ylim(0.008, 1.8)
# ax.set_xlabel(r'Peak flux (erg s$^{-1}$ cm$^{-2}$)', fontsize=12)
# ax.set_ylabel('CCDF', fontsize=12)
# ax.set_title('CCDF of Mg\u202fII h pre-flare transients', fontsize=12)
# ax.legend(fontsize=9.5, frameon=True, framealpha=0.9, loc='upper right')
# ax.tick_params(which='both', direction='in', top=True, right=True)
# ax.grid(True, which='both', alpha=0.25, lw=0.5)
# ax.annotate(f'$N={N}$,  $n_1={n1}$,  $n_2={n2}$',
#             xy=(0.97, 0.06), xycoords='axes fraction',
#             ha='right', va='bottom', fontsize=8.5, color='#5F5E5A')

# plt.tight_layout()
# plt.savefig('ccdf_corrected.pdf', dpi=300, bbox_inches='tight')



import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar

# df = pd.read_csv('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Comapre_flares/pre_flare_intensity_distribution/Flux_distrib/Peak_Flux.csv#')
# data = np.sort(df['intensity'].values)
data = np.sort(data[data > 0.36])
EA, EB = 0.36, data.max()
N = len(data)
xb = 1.03611
seg1 = data[data < xb];  n1 = len(seg1)
seg2 = data[data >= xb]; n2 = len(seg2)

def logL(a, x, xlo, xhi):
    if a <= 1: return -np.inf
    norm = (xhi**(1-a) - xlo**(1-a)) / (1-a)
    return -a*np.sum(np.log(x)) - len(x)*np.log(norm)

def mle(x, xlo, xhi):
    return minimize_scalar(lambda a: -logL(a,x,xlo,xhi),
                           bounds=(1.001,8), method='bounded').x

a1 = mle(seg1, EA, xb)
a2 = mle(seg2, xb, EB)

np.random.seed(42)
b1, b2 = [], []
for _ in range(2000):
    b1.append(mle(np.random.choice(seg1,n1,replace=True), EA, xb))
    b2.append(mle(np.random.choice(seg2,n2,replace=True), xb, EB))
b1, b2 = np.array(b1), np.array(b2)
a1_lo, a1_hi = np.percentile(b1,16), np.percentile(b1,84)
a2_lo, a2_hi = np.percentile(b2,16), np.percentile(b2,84)

# ── Straight-line (unbounded) PL for visual representation ──
# Seg1 anchored at (EA, n1/N): CCDF = (n1/N) * (x/EA)^{-(a1-1)}
# Seg2 anchored at (xb, n2/N): CCDF = (n2/N) * (x/xb)^{-(a2-1)}
x1 = np.logspace(np.log10(EA), np.log10(xb), 400)
x2 = np.logspace(np.log10(xb), np.log10(EB), 400)

c1 = (n1/N) * (x1/EA)**(-(a1-1))
c2 = (n2/N) * (x2/xb)**(-(a2-1))

lo1 = np.percentile([(n1/N)*(x1/EA)**(-(b-1)) for b in b1], 16, axis=0)
hi1 = np.percentile([(n1/N)*(x1/EA)**(-(b-1)) for b in b1], 84, axis=0)
lo2 = np.percentile([(n2/N)*(x2/xb)**(-(b-1)) for b in b2], 16, axis=0)
hi2 = np.percentile([(n2/N)*(x2/xb)**(-(b-1)) for b in b2], 84, axis=0)

# print("Straight-line check:")
# print(f"  Seg1 at EA={EA:.3f}: {c1[0]:.5f}  (expected n1/N={n1/N:.5f})")
# print(f"  Seg1 at xb={xb:.3f}: {c1[-1]:.5f}")
# print(f"  Seg2 at xb={xb:.3f}: {c2[0]:.5f}  (expected n2/N={n2/N:.5f})")
# print(f"  Seg2 at EB={EB:.3f}: {c2[-1]:.5f}")
# print(f"\n  Seg1 is straight line on log-log: True (by construction)")
# print(f"  Seg2 is straight line on log-log: True (by construction)")
# print(f"  Small overlap region near xb is intentional -- shows visual continuity")



x_emp = np.sort(data_)
N=len(data_)
print(N)
ccdf_emp = np.arange(N,0,-1)/N

fig, ax = plt.subplots(figsize=(7, 5))

ax.plot(x_emp, ccdf_emp, '+', color='#2C2C2A', ms=5, mew=0.9,
        zorder=5, label='SUIT Mg II h transients')

ax.fill_between(x1, lo1*1.2, hi1*1.2, color='#185FA5', alpha=0.13, linewidth=0)
ax.fill_between(x2, lo2*0.65, hi2*0.65, color='#D85A30', alpha=0.13, linewidth=0)

ax.plot(x1, c1*1.2, color='#185FA5', lw=2.0, zorder=4, label=rf'$\alpha_1 = {a1:.2f}^{{+{a1_hi-a1:.2f}}}_{{-{a1-a1_lo:.2f}}}$')
ax.plot(x2, c2*0.65, color='#D85A30', lw=2.0, zorder=4, label=rf'$\alpha_2 = {a2:.2f}^{{+{a2_hi-a2:.2f}}}_{{-{a2-a2_lo:.2f}}}$')

ax.axvline(EA, color="#1725BA", lw=1.4, ls='--', zorder=3,
           label=rf'$x_a = {EA:.3f}$ erg s$^{{-1}}$ cm$^{{-2}}$')

ax.axvline(xb, color='#BA7517', lw=1.4, ls='--', zorder=3,
           label=rf'$x_b = {xb:.3f}$ erg s$^{{-1}}$ cm$^{{-2}}$')

ax.set_xscale('log'); ax.set_yscale('log')
# ax.set_xlim(0.28, 4.0); ax.set_ylim(0.008, 1.8)
ax.set_xlabel(r'Peak flux (erg s$^{-1}$ cm$^{-2}$)', fontsize=12)
ax.set_ylabel('CCDF', fontsize=12)
ax.set_title('CCDF of Mg II h pre-flare transients', fontsize=12)
ax.legend(fontsize=10, frameon=True, framealpha=0.9,edgecolor='#D3D1C7', loc='lower left')
ax.tick_params(which='both', direction='in', top=True, right=True)
ax.grid(True, which='both', alpha=0.25, lw=0.5)
# ax.annotate(rf'$N={N}$,  $n_1={n1}$,  $n_2={n2}$', xy=(0.59, 0.06), xycoords='axes fraction', ha='right', va='bottom', fontsize=8.5, color='#5F5E5A')

plt.tight_layout()
plt.savefig('ccdf_straight_lines.pdf', dpi=300, bbox_inches='tight')
plt.savefig('ccdf_straight_lines.png', dpi=300, bbox_inches='tight')
print("Saved.")


