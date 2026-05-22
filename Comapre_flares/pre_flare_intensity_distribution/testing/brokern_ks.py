import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

data = np.loadtxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Comapre_flares/pre_flare_intensity_distribution/Flux_distrib/Peak_Flux.csv')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

E = np.sort(data)
N = len(E)
ccdf_all = np.arange(N, 0, -1) / N

def mle_alpha(x, xmin):
    x = x[x >= xmin]
    n = len(x)
    return 1.0 + n / np.sum(np.log(x / xmin))

def ks_distance_fixed_xmin(x, xmin):
    x = x[x >= xmin]
    n = len(x)
    if n < 15:
        return np.nan, np.nan
    alpha    = mle_alpha(x, xmin)
    ccdf_emp = np.arange(n, 0, -1) / n
    ccdf_mod = (x / xmin) ** (-(alpha - 1))
    D        = np.max(np.abs(ccdf_emp - ccdf_mod))
    return D, alpha

# ── Step 1: Fix Emin from physics (bottom of real power-law region) ───────────
# From your CCDF plot the data starts declining consistently around 0.05
# The KS scan WITHIN this fixed lower segment finds the best α1
Emin_fixed = 0.05     # adjust: where the flat plateau ends

# Scan for best α1 only within [Emin_fixed, some upper cap]
# Upper cap keeps the break region out of this fit
upper_cap  = 0.60     # scan seg1 only below the visible bend

seg1_candidates = np.unique(E[(E >= Emin_fixed) & (E < upper_cap)])
ks1, a1s, emins1 = [], [], []

for emin in seg1_candidates:
    seg = E[(E >= emin) & (E < upper_cap)]
    if len(seg) < 15:
        continue
    D, a = ks_distance_fixed_xmin(seg, emin)
    if np.isnan(D):
        continue
    ks1.append(D); a1s.append(a); emins1.append(emin)

ks1    = np.array(ks1)
a1s    = np.array(a1s)
emins1 = np.array(emins1)

idx1      = np.argmin(ks1)
Emin_opt  = emins1[idx1]
alpha1    = a1s[idx1]
sigma1    = (alpha1 - 1) / np.sqrt(np.sum((E >= Emin_opt) & (E < upper_cap)))

print(f"Seg 1 — Emin  = {Emin_opt:.4f}")
print(f"Seg 1 — α     = {alpha1:.3f} ± {sigma1:.3f}")
print(f"Seg 1 — N     = {np.sum((E >= Emin_opt) & (E < upper_cap))}")

# ── Step 2: Scan for break energy above seg1 ──────────────────────────────────
# Search only above upper_cap so it doesn't overlap with seg1
break_candidates = np.unique(E[(E >= upper_cap) & (E < E[-15])])
ks2, a2s, ebreaks = [], [], []

for eb in break_candidates:
    seg = E[E >= eb]
    if len(seg) < 15:
        continue
    D, a = ks_distance_fixed_xmin(seg, eb)
    if np.isnan(D):
        continue
    ks2.append(D); a2s.append(a); ebreaks.append(eb)

ks2     = np.array(ks2)
a2s     = np.array(a2s)
ebreaks = np.array(ebreaks)

idx2    = np.argmin(ks2)
Ebreak  = ebreaks[idx2]
alpha2  = a2s[idx2]
sigma2  = (alpha2 - 1) / np.sqrt(np.sum(E >= Ebreak))

print(f"\nSeg 2 — Ebreak = {Ebreak:.4f}")
print(f"Seg 2 — α      = {alpha2:.3f} ± {sigma2:.3f}")
print(f"Seg 2 — N      = {np.sum(E >= Ebreak)}")

# ── Step 3: Bootstrap ─────────────────────────────────────────────────────────
def bootstrap_alpha(seg, xmin, n_boot=2000):
    rng = np.random.default_rng(42)
    out = []
    for _ in range(n_boot):
        s = np.sort(rng.choice(seg, len(seg), replace=True))
        s = s[s >= xmin]
        if len(s) < 5: continue
        out.append(mle_alpha(s, xmin))
    return np.mean(out), np.std(out)

seg1 = E[(E >= Emin_opt) & (E < Ebreak)]
seg2 = E[E >= Ebreak]

a1m, a1s_boot = bootstrap_alpha(seg1, Emin_opt)
a2m, a2s_boot = bootstrap_alpha(seg2, Ebreak)

print(f"\nBootstrap: α1 = {a1m:.3f} ± {a1s_boot:.3f}")
print(f"Bootstrap: α2 = {a2m:.3f} ± {a2s_boot:.3f}")

# ── Step 4: KS diagnostic plots ───────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(11, 4))

axes[0].plot(emins1, ks1, color="black", lw=1.2)
axes[0].axvline(Emin_opt, color="#2E86AB", ls="--", lw=1.5,
                label=f"Emin = {Emin_opt:.3f}")
axes[0].set_xscale("log")
axes[0].set_xlabel(r"$E_{\min}$ candidate", fontsize=11)
axes[0].set_ylabel("KS distance D", fontsize=11)
axes[0].set_title("Seg 1 — KS scan (below break)", fontsize=11)
axes[0].legend(fontsize=9)
axes[0].grid(True, which="both", ls=":", lw=0.4, alpha=0.5)

axes[1].plot(ebreaks, ks2, color="black", lw=1.2)
axes[1].axvline(Ebreak, color="#E07B39", ls="--", lw=1.5,
                label=f"Ebreak = {Ebreak:.3f}")
axes[1].set_xscale("log")
axes[1].set_xlabel(r"$E_{\rm break}$ candidate", fontsize=11)
axes[1].set_ylabel("KS distance D", fontsize=11)
axes[1].set_title("Seg 2 — KS scan (above break)", fontsize=11)
axes[1].legend(fontsize=9)
axes[1].grid(True, which="both", ls=":", lw=0.4, alpha=0.5)

plt.tight_layout()
plt.savefig("ks_scans_v2.pdf", dpi=180, bbox_inches="tight")
plt.show()