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
    x = np.array(data['Img_count'], dtype=float)#/area
    
    case_id=os.path.basename(f)[0:3]
    flt=os.path.basename(f)[-8:-4]
    #print(flt)
    
    if case_id=='c04':
        t_start =  np.datetime64("2024-07-10T03:59:00")
        t_end   =  np.datetime64("2024-07-10T05:44:00")
        norm_v=norm[0]

    elif case_id=='c05':
        t_start =  np.datetime64("2024-07-10T13:37:00")
        t_end   =  np.datetime64("2024-07-10T15:25:00")
        norm_v=norm[1]

    elif case_id=='c06':
        t_start =  np.datetime64("2024-10-08T23:56:00")
        t_end   =  np.datetime64("2024-10-09T01:25:00")
        norm_v=norm[2]

    elif case_id=='c07':
        t_start =  np.datetime64("2024-11-01T00:16:00")
        t_end   =  np.datetime64("2024-11-01T02:05:00")
        norm_v=norm[3]

    elif case_id=='c08':
        t_start =  np.datetime64("2024-11-01T12:31:00")
        t_end   =  np.datetime64("2024-11-01T14:18:00")
        norm_v=norm[4]

    elif case_id=='c09':
        t_start =  np.datetime64("2024-11-12T22:22:00")
        t_end   =  np.datetime64("2024-11-13T00:10:00")
        norm_v=norm[5]

    elif case_id=='c10':
        t_start =  np.datetime64("2024-11-13T15:08:00")
        t_end   =  np.datetime64("2024-11-13T16:57:00")
        norm_v=norm[6]
    
    mask = (t >= t_start) & (t <= t_end)
    t_cut = t[mask]
    img_count_cut  = x[mask]/norm_v
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

    # for i in range(len(peaks)):
    #     plt.axvline(t_cut[peaks[i]], color='b',alpha=0.2)
    for k in range(len(h_pks)):
        plt.axvline(h_pks[k], color='r',alpha=0.2)
        # print(h_pks[k],h_counts[k])
    plt.scatter(h_pks,h_counts*1e7,marker='+',c='magenta')
    #plt.yscale('log')
    plt.title(f'{case_id}')
    plt.axhline(mid_cut,color='r',linestyle='-',label='Median')
    plt.axhline(threshold,color='g', linestyle=':',label='Threshold')
    plt.legend()
    plt.savefig(f'pks/{case_id}_{flt}_local_peaks.png',dpi=300)
    plt.close()
    # print('sigma_thresh : ', np.std(img_count_cut)*0.5)


# Convert to datetime64 if not already
t_peaks = np.array(t_peaks, dtype='datetime64[ms]')

# tolerance
tol = np.timedelta64(30, 's')

matched_ints = []
all_c = np.concatenate(all_ints)
all_t = np.concatenate(all_times)

print('number peaks: ',len(all_c),len(t_peaks))

# for dt_arr, pk_arr in zip(all_t, all_c):
#     dt_arr = np.array(dt_arr, dtype='datetime64[ms]')
    
#     for t_x in t_peaks:
#         mask = np.abs(dt_arr - t_x) <= tol
#         matched_ints.extend(pk_arr[mask])

# matched_ints = np.array(matched_ints)
#------------------------------------
# for dt_arr, pk_arr in zip(all_times, all_ints):
#     dt_arr = np.array(dt_arr, dtype='datetime64[ms]')
    
#     for t_x in t_peaks:
#         mask = np.abs(dt_arr - t_x) <= tol
#         matched_ints.extend(pk_arr[mask])

# matched_ints = np.array(matched_ints)


E = np.array(all_c)


E = E[E > 0]
E = np.sort(E)

# print('---------->',len(E))
# # E = E[60:]
# def mle_alpha(E, Emin):
#     return 1.0 + len(E) / np.sum(np.log(E / Emin))

# def empirical_ccdf(E):
#     n = len(E)
#     return np.arange(n, 0, -1) / n

# def model_ccdf(E, Emin, alpha):
#     return (E / Emin) ** (-(alpha - 1))

# def ks_distance(E, Emin):
#     Efit = E[E >= Emin]
#     n = len(Efit)

#     # Require enough statistics
#     if n < 30:
#         return np.nan, np.nan

#     alpha = mle_alpha(Efit, Emin)

#     ccdf_emp = empirical_ccdf(Efit)
#     ccdf_mod = model_ccdf(Efit, Emin, alpha)

#     D = np.max(np.abs(ccdf_emp - ccdf_mod))
#     return D, alpha
# Emins = np.unique(E)
# ks_vals = []
# alphas = []

# for Emin in Emins:
#     D, a = ks_distance(E, Emin)
#     ks_vals.append(D)
#     alphas.append(a)


# E_sorted=E
# Eb1=1.3e7
# Eb2=4.5e7

# from scipy.optimize import minimize_scalar

# def mle_truncated(E, Emin, Emax):
#     Eseg = E[(E >= Emin) & (E < Emax)]
#     N = len(Eseg)
    
#     if N < 20:
#         return np.nan
    
#     sum_logE = np.sum(np.log(Eseg))
    
#     def neg_logL(alpha):
#         if alpha <= 1:
#             return np.inf
        
#         C = (alpha - 1) / (Emin**(1-alpha) - Emax**(1-alpha))
#         logL = N*np.log(C) - alpha*sum_logE
#         return -logL  # minimize
    
#     res = minimize_scalar(neg_logL, bounds=(1.01, 5), method='bounded')
    
#     return res.x

# a1 = mle_truncated(E, E.min(), Eb1)
# a2 = mle_truncated(E, Eb1, Eb2)
# a3 = mle_truncated(E, Eb2, E.max())

# print(a1, a2, a3)

# def ccdf_segment(E, Emin, alpha):
#     return (E / Emin)**(-(alpha - 1))



# plt.loglog(E_sorted, empirical_ccdf(E_sorted), '.', alpha=0.5)

# m1 = E_sorted < Eb1
# m2 = (E_sorted >= Eb1) & (E_sorted < Eb2)
# m3 = E_sorted >= Eb2

# plt.loglog(E_sorted[m1],
#            ccdf_segment(E_sorted[m1], E_sorted[m1][0], a1),
#            '-', label=f'α1={a1:.2f}')

# plt.loglog(E_sorted[m2],
#            ccdf_segment(E_sorted[m2], Eb1, a2),
#            '-', label=f'α2={a2:.2f}')

# plt.loglog(E_sorted[m3],
#            ccdf_segment(E_sorted[m3], Eb2, a3),
#            '-', label=f'α3={a3:.2f}')
# plt.legend()
# plt.show()
#---------------------------
from scipy.optimize import minimize_scalar
import numpy as np

def mle_truncated(E, Emin, Emax):
    Eseg = E[(E >= Emin) & (E < Emax)]
    N = len(Eseg)
    
    if N < 20:
        return np.nan, -np.inf
    
    sum_logE = np.sum(np.log(Eseg))
    
    def neg_logL(alpha):
        if alpha <= 1:
            return np.inf
        
        C = (alpha - 1) / (Emin**(1-alpha) - Emax**(1-alpha))
        logL = N*np.log(C) - alpha*sum_logE
        return -logL
    
    res = minimize_scalar(neg_logL, bounds=(1.01, 5), method='bounded')
    
    return res.x, -res.fun  # alpha, logL

def total_logL(E, Eb1, Eb2):
    Emin = E.min()
    Emax = E.max()
    
    a1, L1 = mle_truncated(E, Emin, Eb1)
    a2, L2 = mle_truncated(E, Eb1, Eb2)
    a3, L3 = mle_truncated(E, Eb2, Emax)
    
    if np.isnan(a1) or np.isnan(a2) or np.isnan(a3):
        return -np.inf, None
    
    return (L1 + L2 + L3), (a1, a2, a3)

E = np.sort(E)

Eb_candidates = np.logspace(np.log10(E.min()), np.log10(E.max()), 30)

best_L = -np.inf
best_params = None

for i in range(5, len(Eb_candidates)-10):
    for j in range(i+5, len(Eb_candidates)-5):
        
        Eb1 = Eb_candidates[i]
        Eb2 = Eb_candidates[j]
        
        L, alphas = total_logL(E, Eb1, Eb2)
        
        if L > best_L:
            best_L = L
            best_params = (Eb1, Eb2, *alphas)

Eb1_opt, Eb2_opt, a1_opt, a2_opt, a3_opt = best_params

print("Eb1, Eb2:", Eb1_opt, Eb2_opt)
print("alphas:", a1_opt, a2_opt, a3_opt)

# a1_opt=1.1
# Eb1=10361258.160823414
# a2_opt=1.6
# Eb2=35368155
# a3_opt=2.8
# print("Eb1, Eb2:", Eb1_opt, Eb2_opt)
# print("alphas:", a1_opt, a2_opt, a3_opt)

E_sorted = np.sort(E)
ccdf = np.arange(len(E_sorted), 0, -1) / len(E_sorted)

plt.loglog(E_sorted, ccdf, '.', alpha=0.4)

# masks
m1 = E_sorted < Eb1_opt
m2 = (E_sorted >= Eb1_opt) & (E_sorted < Eb2_opt)
m3 = E_sorted >= Eb2_opt

# # index of Eb1 and Eb2 in sorted array
# idx1 = np.searchsorted(E_sorted, Eb1, side='left')
# idx2 = np.searchsorted(E_sorted, Eb2, side='left')

# # anchor CCDF values
# scale1 = ccdf[0]        # first segment anchored at start (~1)
# scale2 = ccdf[idx1]     # anchor at Eb1
# scale3 = ccdf[idx2]     # anchor at Eb2

# # Segment 1
# E1 = E_sorted[m1]
# ccdf1 = scale1 * (E1 / E1[0])**(-(a1_opt - 1))

# # Segment 2
# E2 = E_sorted[m2]
# ccdf2 = scale2 * (E2 / Eb1)**(-(a2_opt - 1))

# # Segment 3
# E3 = E_sorted[m3]
# ccdf3 = scale3 * (E3 / Eb2)**(-(a3_opt - 1))



E1 = np.logspace(np.log10(E_sorted.min()), np.log10(Eb1), 100)
E2 = np.logspace(np.log10(Eb1), np.log10(Eb2), 100)
E3 = np.logspace(np.log10(Eb2), np.log10(E_sorted.max()), 100)

idx1 = np.searchsorted(E_sorted, Eb1)
idx2 = np.searchsorted(E_sorted, Eb2)

scale1 = ccdf[0]
# Segment 1
E1 = np.logspace(np.log10(E_sorted.min()), np.log10(Eb1), 100)
ccdf1 = (E1 / E1[0])**(-(a1_opt - 1))

# Value at Eb1 from segment 1
ccdf_Eb1 = ccdf1[-1]

# Segment 2 (continuous)
E2 = np.logspace(np.log10(Eb1), np.log10(Eb2), 100)
ccdf2 = ccdf_Eb1 * (E2 / Eb1)**(-(a2_opt - 1))

# Value at Eb2 from segment 2
ccdf_Eb2 = ccdf2[-1]

# Segment 3 (continuous)
E3 = np.logspace(np.log10(Eb2), np.log10(E_sorted.max()), 100)
ccdf3 = ccdf_Eb2 * (E3 / Eb2)**(-(a3_opt - 1))

plt.figure(figsize=(6,5))

plt.loglog(E_sorted, ccdf, '.', alpha=0.4)

plt.loglog(E1, ccdf1, '-', lw=2, label=f'α1={a1_opt:.2f}')
plt.loglog(E2, ccdf2, '-', lw=2, label=f'α2={a2_opt:.2f}')
plt.loglog(E3, ccdf3, '-', lw=2, label=f'α3={a3_opt:.2f}')

plt.axvline(Eb1, color='k', linestyle='--')
plt.axvline(Eb2, color='k', linestyle='--')

plt.legend()
plt.savefig('pow_fit.png',dpi=300)
plt.show()


# plt.figure(figsize=(6,5))

# # data
# plt.loglog(E_sorted, ccdf, '.', alpha=0.4)

# # model segments
# plt.loglog(E1, ccdf1, '-', label=f'α1={a1_opt:.2f}')
# plt.loglog(E2, ccdf2, '-', label=f'α2={a2_opt:.2f}')
# plt.loglog(E3, ccdf3, '-', label=f'α3={a3_opt:.2f}')

# # breakpoints
# plt.axvline(Eb1, color='k', linestyle='--')
# plt.axvline(Eb2, color='k', linestyle='--')

# plt.xlabel('E')
# plt.ylabel('CCDF')
# plt.legend()
# plt.grid(True, which='both', alpha=0.3)
# plt.savefig('pow_fit.png',dpi=300)
# plt.show()

# def ccdf_model(E, Emin, alpha):
#     return (E / Emin)**(-(alpha - 1))

# plt.loglog(E_sorted[m1],
#            ccdf_model(E_sorted[m1], E_sorted[m1][0], a1_opt),
#            '-', label=f'α1={a1_opt:.2f}')

# idx = np.searchsorted(E_sorted, Eb1)
# scale = ccdf[idx]

# ccdf_model = scale * (E_segment / Eb1)**(-(alpha2 - 1))

# plt.loglog(E_sorted[m2],
#            ccdf_model(E_sorted[m2], Eb1_opt, a2_opt),
#            '-', label=f'α2={a2_opt:.2f}')

# plt.loglog(E_sorted[m3],
#            ccdf_model(E_sorted[m3], Eb2_opt, a3_opt),
#            '-', label=f'α3={a3_opt:.2f}')

# plt.axvline(Eb1_opt, color='k', linestyle='--')
# plt.axvline(Eb2_opt, color='k', linestyle='--')

# plt.legend()
# plt.xlabel('E')
# plt.ylabel('CCDF')
# plt.savefig('ccd_pow_law.png',dpi=300)
# plt.show()


# def mle_truncated(E, Emin, Emax):
#     Eseg = E[(E >= Emin) & (E < Emax)]
#     N = len(Eseg)
    
#     if N < 20:
#         return np.nan, -np.inf
    
#     sum_logE = np.sum(np.log(Eseg))
    
#     def neg_logL(alpha):
#         if alpha <= 1:
#             return np.inf
        
#         C = (alpha - 1) / (Emin**(1-alpha) - Emax**(1-alpha))
#         logL = N*np.log(C) - alpha*sum_logE
#         return -logL
    
#     res = minimize_scalar(neg_logL, bounds=(1.01, 5), method='bounded')
    
#     return res.x, -res.fun  # alpha, logL

# def total_logL(E, Eb2):
#     Emin = E.min()
#     Emax = E.max()
    
#     a1, L1 = mle_truncated(E, Emin, Eb2)
#     a3, L3 = mle_truncated(E, Eb2, Emax)

    
#     return (L1  + L3), (a1, a3)

# E = np.sort(E)

# Eb_candidates = np.logspace(np.log10(E.min()), np.log10(E.max()), 30)

# best_L = -np.inf
# best_params = None

# for j in range(5, len(Eb_candidates)-10):
#     # for j in range(i+5, len(Eb_candidates)-5):
        
#     Eb2 = Eb_candidates[j]
    
#     L, alphas = total_logL(E, Eb2)
    
#     if L > best_L:
#         best_L = L
#         best_params = (Eb2, *alphas)

# Eb2_opt, a1_opt, a3_opt = best_params

# print(" Eb2:", Eb2_opt)
# print("alphas:", a1_opt ,a3_opt)


# E_sorted = np.sort(E)
# ccdf = np.arange(len(E_sorted), 0, -1) / len(E_sorted)

# plt.loglog(E_sorted, ccdf, '.', alpha=0.4)

# # masks
# m1 = E_sorted < Eb2_opt
# # m2 = (E_sorted >= Eb1_opt) & (E_sorted < Eb2_opt)
# m3 = E_sorted >= Eb2_opt

# def ccdf_model(E, Emin, alpha):
#     return (E / Emin)**(-(alpha - 1))

# plt.loglog(E_sorted[m1],
#            ccdf_model(E_sorted[m1], E_sorted[m1][0], a1_opt),
#            '-', label=f'α1={a1_opt:.2f}')


# plt.loglog(E_sorted[m3],
#            ccdf_model(E_sorted[m3], Eb2_opt, a3_opt),
#            '-', label=f'α3={a3_opt:.2f}')

# # plt.axvline(Eb1_opt, color='k', linestyle='--')
# plt.axvline(Eb2_opt, color='k', linestyle='--')

# plt.legend()
# plt.xlabel('E')
# plt.ylabel('CCDF')
# plt.show()