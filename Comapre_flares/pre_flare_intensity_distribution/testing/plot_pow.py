import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Comapre_flares/pre_flare_intensity_distribution/Flux_distrib/Peak_Flux.csv')
data=np.sort(data)
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

# 68% bands (recompute from bootstrap arrays b1, b2)
lo1 = np.percentile([ccdf_trunc(x1,b,EA,xb,n1/N) for b in b1], 16, axis=0)
hi1 = np.percentile([ccdf_trunc(x1,b,EA,xb,n1/N) for b in b1], 84, axis=0)
lo2 = np.percentile([ccdf_trunc(x2,b,xb,EB,n2/N) for b in b2], 16, axis=0)
hi2 = np.percentile([ccdf_trunc(x2,b,xb,EB,n2/N) for b in b2], 84, axis=0)

# ── Plot ──────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 5))

x_emp    = np.sort(data)
ccdf_emp = np.arange(N, 0, -1) / N
ax.plot(x_emp, ccdf_emp, '+', color='#2C2C2A', ms=5, mew=0.9,
        zorder=5, label='SUIT Mg\u202fII h transients')

ax.fill_between(x1, lo1, hi1, color='#185FA5', alpha=0.15, linewidth=0)
ax.fill_between(x2, lo2, hi2, color='#D85A30', alpha=0.15, linewidth=0)

ax.plot(x1, c1, color='#185FA5', lw=2.0, zorder=4,
        label=rf'$\alpha_1={a1:.2f}_{{-{a1-a1_lo:.2f}}}^{{+{a1_hi-a1:.2f}}}$'
              r'  ($x < x_b$)')
ax.plot(x2, c2, color='#D85A30', lw=2.0, zorder=4,
        label=rf'$\alpha_2={a2:.2f}_{{-{a2-a2_lo:.2f}}}^{{+{a2_hi-a2:.2f}}}$'
              r'  ($x \geq x_b$)')

ax.axvline(xb, color='#BA7517', lw=1.4, ls='--', zorder=3,
           label=rf'$x_b = {xb:.3f}$ erg s$^{{-1}}$ cm$^{{-2}}$')

ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(0.28, 4.0); ax.set_ylim(0.008, 1.8)
ax.set_xlabel(r'Peak flux (erg s$^{-1}$ cm$^{-2}$)', fontsize=12)
ax.set_ylabel('CCDF', fontsize=12)
ax.set_title('CCDF of Mg\u202fII h pre-flare transients', fontsize=12)
ax.legend(fontsize=9.5, frameon=True, framealpha=0.9, loc='upper right')
ax.tick_params(which='both', direction='in', top=True, right=True)
ax.grid(True, which='both', alpha=0.25, lw=0.5)
ax.annotate(f'$N={N}$,  $n_1={n1}$,  $n_2={n2}$',
            xy=(0.97, 0.06), xycoords='axes fraction',
            ha='right', va='bottom', fontsize=8.5, color='#5F5E5A')

plt.tight_layout()
plt.savefig('ccdf_corrected.pdf', dpi=300, bbox_inches='tight')