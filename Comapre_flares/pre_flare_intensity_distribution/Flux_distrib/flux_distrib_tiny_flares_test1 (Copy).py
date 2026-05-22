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
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
# set_pub_style()
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution, curve_fit
from scipy.stats import kstest
from scipy.optimize import minimize_scalar
from itertools import product


scol =sns.color_palette("colorblind")

files = sorted(glob.glob("csv_files/*NB04.csv"))
time_fls=sorted(glob.glob("x-ray_peaks/*.csv"))
norm =np.genfromtxt('norm_vals.csv',delimiter=',')


all_times = []
all_ints = []
t_peaks=[]
for t in time_fls:
    tp=np.genfromtxt(t, delimiter=',', names=True, dtype=None, encoding='utf-8')
    t_peaks.append(tp['date_time'])

t_peaks=np.concatenate(t_peaks)
# print(t_peaks)
l=0
for f in files:

    data = np.genfromtxt(f, delimiter=',', names=True, dtype=None, encoding='utf-8')
    hpks=(np.genfromtxt(time_fls[l], delimiter=',', names=True, dtype=None, encoding='utf-8'))['date_time']
    h_counts=(np.genfromtxt(time_fls[l], delimiter=',', names=True, dtype=None, encoding='utf-8'))['helios_counts']
    l+=1
    h_pks=np.array(hpks,dtype='datetime64[ms]')
    # Convert time to datetime64
    t = np.array(data['Date'], dtype='datetime64[ms]')
    area=np.array(data['area'],dtype=float)
    # x = np.array(data['diff_count'], dtype=float)#/area
    x = np.array(data['Img_count'], dtype=float)#*(area*11.444e-6)#(np.pi/3600/180)**2*0.698*0.698)
    # print(x)
    case_id=os.path.basename(f)[0:3]
    flt=os.path.basename(f)[-8:-4]
    #print(flt)
    
    if case_id=='c04':
        t_start =  np.datetime64("2024-07-10T03:59:00")
        t_end   =  np.datetime64("2024-07-10T05:44:00")
        m=1.87e3
        c=3.39e6
        

    elif case_id=='c05':
        t_start =  np.datetime64("2024-07-10T13:37:00")
        t_end   =  np.datetime64("2024-07-10T15:25:00")
        m=1.87e3
        c=3.39e6

    elif case_id=='c06':
        t_start =  np.datetime64("2024-10-08T23:56:00")
        t_end   =  np.datetime64("2024-10-09T01:25:00")
        m=2.32e3
        c=-2.77e6
        

    elif case_id=='c07':
        t_start =  np.datetime64("2024-11-01T00:16:00")
        t_end   =  np.datetime64("2024-11-01T02:05:00")
        m=2.84e3
        c=-4.78e6
        

    elif case_id=='c08':
        t_start =  np.datetime64("2024-11-01T12:31:00")
        t_end   =  np.datetime64("2024-11-01T14:18:00")
        m=2.84e3
        c=-4.78e6

    elif case_id=='c09':
        t_start =  np.datetime64("2024-11-12T22:22:00")
        t_end   =  np.datetime64("2024-11-13T00:10:00")
        m=3.28e3
        c=-4.43e6

    elif case_id=='c10':
        t_start =  np.datetime64("2024-11-13T15:08:00")
        t_end   =  np.datetime64("2024-11-13T16:57:00")
        m=3.28e3
        c=-4.43e6
    
    mask = (t >= t_start) & (t <= t_end)
    t_cut = t[mask]
    img_count_cut  = (x[mask]*m+c*area[mask])*((np.pi/3600/180)**2*0.698*0.698)
    # print(area[mask])
    mid_cut=np.median(img_count_cut)

    fig,ax=plt.subplots(figsize=(14,6))

    threshold=0#2.5e6 #
    peaks_=find_peaks(img_count_cut, height=threshold)
    dt_=np.array(t_cut[peaks_[0]])
    pk_=np.array(img_count_cut[peaks_[0]])
    for i in range(len(peaks_[0])):
        plt.axvline(t_cut[peaks_[0][i]], color='g',alpha=0.2)

    peaks,props = find_peaks(img_count_cut, prominence=0)#,width=3)
    dt=dt_#np.array(t_cut[peaks])
    pk=pk_#np.array(img_count_cut[peaks])

    all_times.append(dt)
    all_ints.append(pk)

    df = pd.DataFrame({"Peak time": dt,  "Peak total count": pk})
    df.to_csv(f"pks/{case_id}_{flt}_peaks.csv", index=False)
    ax.plot(t_cut,img_count_cut)

    for k in range(len(h_pks)):
        plt.axvline(h_pks[k], color='r',alpha=0.2)
        # print(h_pks[k],h_counts[k])
    plt.scatter(h_pks,h_counts,marker='+',c='magenta')
    #plt.yscale('log')
    plt.title(f'{case_id}')
    plt.axhline(mid_cut,color='r',linestyle='-',label='Median')
    plt.axhline(threshold,color='g', linestyle=':',label='Threshold')
    plt.legend()
    plt.yscale('log')
    plt.savefig(f'pks/{case_id}_{flt}_local_peaks.png',dpi=300)
    plt.close()


# Convert to datetime64 if not already
t_peaks = np.array(t_peaks, dtype='datetime64[ms]')
tol = np.timedelta64(30, 's')
matched_ints = []
all_c = np.concatenate(all_ints)
all_t = np.concatenate(all_times)

print('number peaks: ',len(all_c),len(t_peaks))

for dt_arr, pk_arr in zip(all_times, all_ints):
    dt_arr = np.array(dt_arr, dtype='datetime64[ms]')
    
    for t_x in t_peaks:
        mask = np.abs(dt_arr - t_x) <= tol
        matched_ints.extend(pk_arr[mask])

matched_ints = np.array(matched_ints)
matched_ints = matched_ints#[matched_ints >= Emin_opt]
matched_ints = np.sort(matched_ints)

E = np.array(all_c)
E = E[E > 0]
E = np.sort(E)


data = E

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
            (-1, 0),   # a1 flat
            (-2,  0),   # a2 moderate
            (-4, 0),   # a3 steep
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

# bic1, aic1, aicc1 = bic_aic(rss_1, 2, n)
# bic2, aic2, aicc2 = bic_aic(rss_2, 4, n)
bic3, aic3, aicc3 = bic_aic(rss_3, 6, n)

print("\n╔══════════════════════════════════════════════════╗")
print("║          MODEL SELECTION SUMMARY                ║")
print("╠══════════════════════════════════════════════════╣")
print(f"║  Model     k    BIC       AIC      AICc        ║")
# print(f"║  1-segment 2  {bic1:8.2f}  {aic1:8.2f}  {aicc1:8.2f}      ║")
# print(f"║  2-segment 4  {bic2:8.2f}  {aic2:8.2f}  {aicc2:8.2f}      ║")
print(f"║  3-segment 6  {bic3:8.2f}  {aic3:8.2f}  {aicc3:8.2f}      ║")
print("╚══════════════════════════════════════════════════╝")

bics   = [bic3] #[bic1, bic2, bic3]
best_k = np.argmin(bics) + 1
# delta  = [b - min(bics) for b in bics]
# print(f"\n  ✔ Best model: {best_k} segment(s)")
# print(f"  ΔBIC vs best:  1-seg={delta[0]:.1f},  "
#       f"2-seg={delta[1]:.1f},  3-seg={delta[2]:.1f}")
# print("  (ΔBIC > 10 = strong evidence against that model)")

# ══════════════════════════════════════════════════════
# 5. PRINT BEST-FIT PARAMETERS
# ══════════════════════════════════════════════════════
lbp1_B, lbp2_B, a1_B, a2_B, a3_B, lC_B = res3.x
# lbp_2,  a1_2,   a2_2, lC_2             = res2.x
# lbp1_B=-2
# lbp2_B=-1.0
# print("\n─── 2-segment parameters ───")
# print(f"  Break energy   : {10**lbp_2:.4f}  erg cm⁻²")
# print(f"  α₁ (below)     : {-a1_2:.3f}")
# print(f"  α₂ (above)     : {-a2_2:.3f}")

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
# ax1.plot(10**lx_fine, 10**y2f, 'b--', lw=1.8,
        #  label=fr'2-seg  $\alpha_1={-a1_2:.2f},\ \alpha_2={-a2_2:.2f}$')

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
# ax1.text(0.03, 0.05,f"Best model (BIC): {best_label}\n",
        #  f"ΔBIC(2-seg)={delta[1]:.1f}  ΔBIC(3-seg)={delta[2]:.1f}",
        #  transform=ax1.transAxes, fontsize=10,
        #  bbox=dict(boxstyle='round', fc='wheat', alpha=0.7))

ax1.set_xscale('log'); ax1.set_yscale('log')
ax1.set_ylabel('CCDF', fontsize=13)
ax1.set_title('CCDF of pre-flare transients — automatic model selection', fontsize=13)
ax1.legend(fontsize=9, loc='lower right')
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