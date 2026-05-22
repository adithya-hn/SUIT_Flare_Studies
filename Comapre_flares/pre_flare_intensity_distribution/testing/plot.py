import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from scipy.integrate import quad


data = np.loadtxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Comapre_flares/pre_flare_intensity_distribution/Flux_distrib/Peak_Flux.csv')
data = np.sort(data[data>0.36])

EA, EB = 0.36, data.max()
N      = len(data)
Ek     = 1.036          # break energy
d1 = data[data <  Ek];  n1 = len(d1)
d2 = data[data >= Ek];  n2 = len(d2)

# ── MLE slopes ───────────────────────────────────────────────────
def mle_seg(x, xlo, xhi):
    r = minimize_scalar(
        lambda a: a*np.sum(np.log(x)) + len(x)*np.log((xhi**(1-a)-xlo**(1-a))/(1-a))
                  if a > 1 else 1e12,
        bounds=(1.001, 8), method='bounded')
    return r.x

a1 = mle_seg(d1, EA, Ek)
a2 = mle_seg(d2, Ek, EB)

# Bootstrap 68% CI
np.random.seed(42)
b1 = np.array([mle_seg(np.random.choice(d1,n1,replace=True),EA,Ek) for _ in range(2000)])
b2 = np.array([mle_seg(np.random.choice(d2,n2,replace=True),Ek,EB) for _ in range(2000)])

# ── Truncated CCDF (the model you actually fitted) ───────────────
def trunc_ccdf(x, a, xlo, xhi, scale):
    norm = (xhi**(1-a) - xlo**(1-a)) / (1-a)
    return scale * np.clip(1 - (x**(1-a) - xlo**(1-a)) / ((1-a)*norm), 0, 1)

x1 = np.logspace(np.log10(EA),  np.log10(Ek),  400)
x2 = np.logspace(np.log10(Ek),  np.log10(EB),  400)
c1_trunc = trunc_ccdf(x1, a1, EA, Ek, n1/N)
c2_trunc = trunc_ccdf(x2, a2, Ek, EB, n2/N)

# Bootstrap bands on truncated CCDF
lo1 = np.percentile([trunc_ccdf(x1,b,EA,Ek,n1/N) for b in b1], 16, axis=0)
hi1 = np.percentile([trunc_ccdf(x1,b,EA,Ek,n1/N) for b in b1], 84, axis=0)
lo2 = np.percentile([trunc_ccdf(x2,b,Ek,EB,n2/N) for b in b2], 16, axis=0)
hi2 = np.percentile([trunc_ccdf(x2,b,Ek,EB,n2/N) for b in b2], 84, axis=0)

# ── Straight-line slope guide (visual aid only) ──────────────────
# Anchored at midpoint of each segment in log-space for clarity
xm1 = np.logspace(np.log10(EA*1.2), np.log10(Ek*0.85), 100)
xm2 = np.logspace(np.log10(Ek*1.1), np.log10(EB*0.9),  100)

# Anchor straight line to pass through the truncated CCDF midpoint
def anchor(x_arr, a, xlo, xhi, n_seg):
    xmid = np.exp((np.log(xlo)+np.log(xhi))/2)
    ccdf_mid = trunc_ccdf(np.array([xmid]), a, xlo, xhi, n_seg/N)[0]
    return ccdf_mid * (x_arr/xmid)**(-(a-1))

sl1 = anchor(xm1, a1, EA, Ek, n1)
sl2 = anchor(xm2, a2, Ek, EB, n2)

# ── Figure ───────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 5))

# Empirical CCDF
x_emp    = np.sort(data)
ccdf_emp = np.arange(N, 0, -1) / N
ax.plot(x_emp, ccdf_emp, '+', color='#2C2C2A', ms=5.5, mew=0.9,
        zorder=6, label=f'SUIT Mg\u202fII h ($N = {N}$)')

# Confidence bands (on truncated model)
ax.fill_between(x1, lo1, hi1, color='#185FA5', alpha=0.15, lw=0)
ax.fill_between(x2, lo2, hi2, color='#D85A30', alpha=0.15, lw=0)

# Truncated CCDF — PRIMARY model curves
da1p=np.percentile(b1,84)-a1; da1m=a1-np.percentile(b1,16)
da2p=np.percentile(b2,84)-a2; da2m=a2-np.percentile(b2,16)
ax.plot(x1, c1_trunc, color='#185FA5', lw=2.2, zorder=4,
        label=rf'$\alpha_1={a1:.2f}^{{+{da1p:.2f}}}_{{-{da1m:.2f}}}$  (truncated MLE)')
ax.plot(x2, c2_trunc, color='#D85A30', lw=2.2, zorder=4,
        label=rf'$\alpha_2={a2:.2f}^{{+{da2p:.2f}}}_{{-{da2m:.2f}}}$  (truncated MLE)')

# Straight lines — VISUAL GUIDE only
ax.plot(xm1, sl1, color='#185FA5', lw=1.2, ls=':', zorder=3, alpha=0.7,
        label='power law slope guide')
ax.plot(xm2, sl2, color='#D85A30', lw=1.2, ls=':', zorder=3, alpha=0.7)

# Break
ax.axvline(Ek, color='#BA7517', lw=1.4, ls='--', zorder=3,
           label=rf'$E_k = {Ek:.3f}$ erg s$^{{-1}}$ cm$^{{-2}}$')

ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(0.27, 4.2); ax.set_ylim(0.006, 1.9)
ax.set_xlabel(r'Peak flux (erg s$^{-1}$ cm$^{-2}$)', fontsize=12)
ax.set_ylabel('CCDF', fontsize=12)
ax.legend(fontsize=9, frameon=True, framealpha=0.93,
          edgecolor='#D3D1C7', loc='upper right')
ax.tick_params(which='both', direction='in', top=True, right=True)
ax.grid(True, which='both', alpha=0.22, lw=0.5)
plt.tight_layout()
plt.savefig('broken_pl_final.pdf', dpi=300, bbox_inches='tight')
plt.show()