"""
STIX Two-Stage Spectral Fitting: Background + Flare (Thermal + Non-Thermal)
============================================================================
Workflow:
  Stage 1 — Fit background interval with thermal model (f_vth) only.
  Stage 2 — Fit each flare time bin with f_vth alone, then f_vth + thick_fn.
             Apply non-thermal onset criterion (delta_ln_L > threshold).

Non-thermal onset criteria (as used in Case 9 / Nov 13 pipeline):
  - Delta log-likelihood > 10–15
  - delta (spectral index) < 7–8
  - e_c (low-energy cutoff) > 5 keV

Author: Adi (adapted for hot onset exploration)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import astropy.units as u
# from sunkit_spex.fitting.objective_functions.optimising_function import minimize_func
from sunkit_spex.fitting.statistics.gaussian import chi_squared
from sunkit_spex.fitting.optimiser_tools.minimiser_tools import scipy_minimize
from sunkit_spex.fitting.io import load_spec
from sunkit_spex.models.physical.thermal import f_vth
from sunkit_spex.models.physical.nonthermal import thick_fn
from sunkit_spex.fitting.fitter import Fitter
from sunkit_spex.fitting.io import load_spec
# ─────────────────────────────────────────────
# 0. USER CONFIGURATION
# ─────────────────────────────────────────────

# STIX science files (replace with your actual paths / remote paths)
SCI_FILES = [
    "stx_*.fits",  # replace with actual file
]


# flare_start = Time("2024-11-01T02:10:20")
# flare_end   = Time("2024-11-01T02:16:00")
# #---------------Input parameters----------------

# spec_file="../../stx_spectrum_2410315184.fits"
# srm_file="../../stx_srm_2410315184.fits"

# case='7_Nov01'


# start_background_time = "2024-11-01T01:59:00"
# end_background_time   = "2024-11-01T02:03:00" 

# Background time interval
BKG_START = "2024-11-01T01:59:00"  # UTC — replace with actual background window
BKG_END   = "2024-11-01T02:03:00"

# Flare time interval (will be split into 20-s bins internally by SunKitSpex)
FLR_START = "2024-11-01T02:10:00"  # replace with flare start
FLR_END   = "2024-11-01T02:20:00"  # replace with flare end

# Energy range for fitting (keV)
E_MIN = 6.0   # keV — below this dominated by low-energy noise / attenuator
E_MAX = 30.0  # keV — adjust based on count statistics

# Non-thermal onset thresholds
DNLNL_THRESH  = 12.0   # delta ln(L) threshold (thermal vs thermal+NT)
DELTA_MAX     = 7.5    # spectral index upper limit for valid NT
EC_MIN_KEV    = 5.0    # minimum valid low-energy cutoff (keV)

# ─────────────────────────────────────────────
# 1. LOAD SPECTRA
# ─────────────────────────────────────────────

spec = load_spec(SCI_FILES)

# Set energy fitting range
spec.energy_fitting_range = [E_MIN, E_MAX] * u.keV

# ─────────────────────────────────────────────
# 2. STAGE 1: FIT BACKGROUND WITH THERMAL ONLY
# ─────────────────────────────────────────────

print("=" * 60)
print("Stage 1: Background thermal fit")
print("=" * 60)

spec.select_time(BKG_START, BKG_END)

# Define thermal model for background
def bkg_model(energy, T=1.0, EM=0.01):
    """
    Thermal bremsstrahlung model for background.
    T  in units of 10^7 K (i.e., T=1.0 → 10 MK)
    EM in units of 10^49 cm^-3
    """
    return f_vth(energy, T, EM)

bkg_result = spec.fit(
    model=bkg_model,
    optimiser=scipy_minimize,
    p0={"T": 1.0, "EM": 0.01},
    bounds={"T": (0.5, 5.0), "EM": (1e-4, 10.0)},
)

T_bkg  = bkg_result.params["T"]
EM_bkg = bkg_result.params["EM"]

print(f"  Background T  = {T_bkg:.3f} (×10 MK)")
print(f"  Background EM = {EM_bkg:.4e} (×10^49 cm^-3)")
print()

# ─────────────────────────────────────────────
# 3. STAGE 2: FLARE BINS — THERMAL ONLY FIRST
# ─────────────────────────────────────────────

print("=" * 60)
print("Stage 2a: Flare bins — thermal-only fit")
print("=" * 60)

spec.select_time(FLR_START, FLR_END)

def flare_thermal_model(energy, T=2.0, EM=0.1):
    return f_vth(energy, T, EM)

th_results = spec.fit_time_bins(
    model=flare_thermal_model,
    optimiser=scipy_minimize,
    p0={"T": 2.0, "EM": 0.1},
    bounds={"T": (0.5, 20.0), "EM": (1e-5, 100.0)},
)

# Log-likelihoods for thermal-only fits
lnL_thermal = np.array([r.lnL for r in th_results])

# ─────────────────────────────────────────────
# 4. STAGE 2b: FLARE BINS — THERMAL + NON-THERMAL
# ─────────────────────────────────────────────

print("=" * 60)
print("Stage 2b: Flare bins — thermal + non-thermal fit")
print("=" * 60)

def flare_th_nth_model(energy, T=2.0, EM=0.1, norm=1e-3, gamma=4.0, e_c=10.0):
    """
    Combined thermal + thick-target non-thermal model.
    T     : temperature in 10^7 K
    EM    : emission measure in 10^49 cm^-3
    norm  : non-thermal normalisation (electrons s^-1 keV^-1 at 50 keV)
    gamma : electron spectral index (delta)
    e_c   : low-energy cutoff (keV)
    """
    th  = f_vth(energy, T, EM)
    nth = thick_fn(energy, norm, gamma, e_c)
    return th + nth

th_nth_results = spec.fit_time_bins(
    model=flare_th_nth_model,
    optimiser=scipy_minimize,
    p0={"T": 2.0, "EM": 0.1, "norm": 1e-3, "gamma": 4.0, "e_c": 10.0},
    bounds={
        "T":     (0.5, 20.0),
        "EM":    (1e-5, 100.0),
        "norm":  (0.0, 1e2),
        "gamma": (2.0, 12.0),
        "e_c":   (5.0, 30.0),
    },
)

lnL_th_nth = np.array([r.lnL for r in th_nth_results])

# ─────────────────────────────────────────────
# 5. NON-THERMAL ONSET CRITERION
# ─────────────────────────────────────────────

print("=" * 60)
print("Stage 3: Non-thermal onset identification")
print("=" * 60)

delta_lnL = lnL_th_nth - lnL_thermal   # should be >= 0

# Extract non-thermal parameters
gamma_arr = np.array([r.params["gamma"] for r in th_nth_results])
ec_arr    = np.array([r.params["e_c"]   for r in th_nth_results])
norm_arr  = np.array([r.params["norm"]  for r in th_nth_results])

# NT onset: all three criteria must be met
nt_onset_mask = (
    (delta_lnL  >= DNLNL_THRESH) &
    (gamma_arr  <= DELTA_MAX)     &
    (ec_arr     >= EC_MIN_KEV)
)

# First bin meeting all criteria
if nt_onset_mask.any():
    nt_onset_bin = np.argmax(nt_onset_mask)
    print(f"  Non-thermal onset at bin index: {nt_onset_bin}")
    print(f"  Onset Δln(L) = {delta_lnL[nt_onset_bin]:.2f}")
    print(f"  Onset δ      = {gamma_arr[nt_onset_bin]:.2f}")
    print(f"  Onset e_c    = {ec_arr[nt_onset_bin]:.2f} keV")
else:
    print("  No non-thermal onset detected with current thresholds.")
    nt_onset_bin = None

# ─────────────────────────────────────────────
# 6. BUILD RESULT ARRAYS
# ─────────────────────────────────────────────

# Mid-times for each bin
bin_times = np.array([
    (r.time_range.start + (r.time_range.end - r.time_range.start) / 2).unix
    for r in th_results
])

# Thermal parameters across bins (use th_nth results for consistency)
T_arr  = np.array([r.params["T"]  for r in th_nth_results])
EM_arr = np.array([r.params["EM"] for r in th_nth_results])

# ─────────────────────────────────────────────
# 7. PUBLICATION-QUALITY FIGURE
# ─────────────────────────────────────────────

fig = plt.figure(figsize=(12, 14))
gs  = GridSpec(4, 1, figure=fig, hspace=0.08, left=0.12, right=0.95,
               top=0.93, bottom=0.07)

ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1], sharex=ax1)
ax3 = fig.add_subplot(gs[2], sharex=ax1)
ax4 = fig.add_subplot(gs[3], sharex=ax1)

t_rel = bin_times - bin_times[0]   # seconds from flare start

# ── Panel 1: Temperature ────────────────────
ax1.plot(t_rel, T_arr * 10, "o-", color="#C0392B", ms=4, lw=1.5, label="T (MK)")
ax1.set_ylabel("Temperature (MK)", fontsize=11)
ax1.legend(fontsize=9, loc="upper left")
ax1.grid(alpha=0.3)

# ── Panel 2: Emission Measure ────────────────
ax2.semilogy(t_rel, EM_arr, "s-", color="#2980B9", ms=4, lw=1.5, label="EM (×10⁴⁹ cm⁻³)")
ax2.set_ylabel("EM (×10⁴⁹ cm⁻³)", fontsize=11)
ax2.legend(fontsize=9, loc="upper left")
ax2.grid(alpha=0.3)

# ── Panel 3: Δln(L) ─────────────────────────
ax3.plot(t_rel, delta_lnL, "D-", color="#27AE60", ms=4, lw=1.5, label="Δln(L)")
ax3.axhline(DNLNL_THRESH, ls="--", color="gray", lw=1.2,
            label=f"Threshold = {DNLNL_THRESH}")
ax3.set_ylabel("Δln(L) [NT − Th]", fontsize=11)
ax3.legend(fontsize=9, loc="upper left")
ax3.grid(alpha=0.3)

# ── Panel 4: NT parameters ──────────────────
ax4_r = ax4.twinx()
ax4.plot(t_rel, gamma_arr, "^-", color="#8E44AD", ms=4, lw=1.5, label="δ (γ)")
ax4_r.plot(t_rel, ec_arr, "v--", color="#E67E22", ms=4, lw=1.5, label="e_c (keV)")
ax4.axhline(DELTA_MAX, ls=":", color="#8E44AD", lw=1.0, alpha=0.7)
ax4_r.axhline(EC_MIN_KEV, ls=":", color="#E67E22", lw=1.0, alpha=0.7)
ax4.set_ylabel("Spectral index δ", fontsize=11)
ax4_r.set_ylabel("Low-energy cutoff e_c (keV)", fontsize=11)
ax4.set_xlabel("Time from flare start (s)", fontsize=11)

lines1, lab1 = ax4.get_legend_handles_labels()
lines2, lab2 = ax4_r.get_legend_handles_labels()
ax4.legend(lines1 + lines2, lab1 + lab2, fontsize=9, loc="upper right")
ax4.grid(alpha=0.3)

# Shade NT onset region in all panels
if nt_onset_bin is not None:
    t_onset = t_rel[nt_onset_bin]
    for ax in [ax1, ax2, ax3, ax4]:
        ax.axvline(t_onset, color="black", ls="-.", lw=1.2, alpha=0.8)
    ax1.text(t_onset + 2, ax1.get_ylim()[1] * 0.95, "NT onset",
             fontsize=8, va="top", ha="left", color="black")

fig.suptitle("STIX Spectral Fitting: Thermal + Non-Thermal (Hot Onset Exploration)",
             fontsize=13, fontweight="bold")

plt.savefig("/mnt/user-data/outputs/stix_two_stage_fit.pdf", dpi=150)
plt.savefig("/mnt/user-data/outputs/stix_two_stage_fit.png", dpi=150)
plt.show()
print("Figures saved.")

# ─────────────────────────────────────────────
# 8. SAVE RESULTS AS NPZ
# ─────────────────────────────────────────────

np.savez(
    "/mnt/user-data/outputs/stix_fit_results.npz",
    bin_times   = bin_times,
    T_arr       = T_arr,
    EM_arr      = EM_arr,
    gamma_arr   = gamma_arr,
    ec_arr      = ec_arr,
    norm_arr    = norm_arr,
    delta_lnL   = delta_lnL,
    nt_onset_mask = nt_onset_mask.astype(int),
    T_bkg       = np.array([T_bkg]),
    EM_bkg      = np.array([EM_bkg]),
)
print("Results saved to stix_fit_results.npz")
