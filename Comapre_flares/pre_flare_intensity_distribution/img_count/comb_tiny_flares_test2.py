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

    for k in range(len(h_pks)):
        plt.axvline(h_pks[k], color='r',alpha=0.2)
        # print(h_pks[k],h_counts[k])
    plt.scatter(h_pks,h_counts*1e7,marker='+',c='magenta')
    #plt.yscale('log')
    plt.title(f'{case_id}')
    plt.axhline(mid_cut,color='r',linestyle='-',label='Median')
    plt.axhline(threshold,color='g', linestyle=':',label='Threshold')
    plt.legend()
    # plt.savefig(f'pks/{case_id}_{flt}_local_peaks.png',dpi=300)
    plt.close()


# Convert to datetime64 if not already
t_peaks = np.array(t_peaks, dtype='datetime64[ms]')
tol = np.timedelta64(30, 's')
matched_ints = []
all_c = np.concatenate(all_ints)
all_t = np.concatenate(all_times)

print('number peaks: ',len(all_c),len(t_peaks))

E = np.array(all_c)
E = E[E > 0]
E = np.sort(E)

# E_scaled = E / E.min()
# E = E_scaled
# Emin = 1

def mle_alpha(Eseg, Emin):
    return 1 + len(Eseg) / np.sum(np.log(Eseg / Emin))

def logL_segment(Eseg, Emin, Emax):
    if len(Eseg) < 20:
        return -np.inf
    
    alpha = mle_alpha(Eseg, Emin)
    
    N = len(Eseg)
    sum_log = np.sum(np.log(Eseg))
    
    C = (alpha - 1) / (Emin**(1-alpha) - Emax**(1-alpha))
    
    return N*np.log(C) - alpha*sum_log
import numpy as np

E = np.sort(E)
E=E[E>1e7]
# candidate breakpoints (log spaced)
Eb_candidates = np.logspace(np.log10(E.min()), np.log10(E.max()), 50)

best_BIC = np.inf
best_params = None

N = len(E)

for i in range(5, len(Eb_candidates)-10):
    
    # Eb1 = Eb_candidates[i]
    Eb2 = Eb_candidates[i]
    
    seg1 = E[E < Eb2]
    # seg2 = E[(E >= Eb1) & (E < Eb2)]
    seg3 = E[E >= Eb2]
    
    if len(seg1) < 20 or len(seg3) < 20:
        continue
    
    L1 = logL_segment(seg1, E.min(), Eb2)
    # L2 = logL_segment(seg2, Eb1, Eb2)
    L3 = logL_segment(seg3, Eb2, E.max())
    
    logL = L1 + L3
    
    k = 5  # α1, α2, α3, Eb1, Eb2
    
    BIC = -2*logL + k*np.log(N)
    
    if BIC < best_BIC:
        best_BIC = BIC
        best_params = ( Eb2)

Eb2_opt = best_params
Eb1_opt=1e7
Eb2_opt=3.4e7
print("Optimal Eb1, Eb2:", Eb1_opt, Eb2_opt)
print("Best BIC:", best_BIC)

# os._exit(1)   
def empirical_ccdf(E):
    n = len(E)
    return np.arange(n, 0, -1) / n

def fit_alphas(E, Eb1, Eb2):
    from scipy.optimize import minimize

    def neg_logL_alpha(params):
        a2, a3 = params
        
        seg1 = E[E < Eb2]
        seg3 = E[E >= Eb2]
        
        def logL(Eseg, Emin, Emax, alpha):
            if alpha <= 1:
                return -np.inf
            N = len(Eseg)
            sum_log = np.sum(np.log(Eseg))
            logC = np.log(alpha - 1) - np.log(Emin**(1-alpha) - Emax**(1-alpha))
            return N*logC - alpha*sum_log
        
        L = (
            logL(seg1, E.b1, Eb2, a2) +
            logL(seg3, Eb2, E.max(), a3)
        )
        
        return -L


    res = minimize(neg_logL_alpha,
                   x0=[1.5,2.9],
                   bounds=[(1,3),(2,4)])   
    return res.x


alphas = []

for _ in range(1000):
    sample = np.random.choice(E, size=len(E), replace=True)
    sample = np.sort(sample)
    
    a2,a3 = fit_alphas(sample, Eb1_opt, Eb2_opt)
    alphas.append([a2,a3])

alphas = np.array(alphas)

a2_mean, a2_std = np.mean(alphas[:,0]), np.std(alphas[:,0])
a3_mean, a3_std = np.mean(alphas[:,1]), np.std(alphas[:,1])
# a3_mean, a3_std = np.mean(alphas[:,2]), np.std(alphas[:,2])

# print(a1_mean, a1_std)
print(a2_mean, a2_std)
print(a3_mean, a3_std)

# a1=a1_mean
a2=a2_mean
a3=a3_mean

E_sorted=E
ccdf = empirical_ccdf(E_sorted)

plt.figure(figsize=(6,5))
plt.loglog(E_sorted, ccdf, '.', label='Data')
# Segment masks
m1 = (E_sorted < Eb1_opt)
m2 = (E_sorted >= Eb1_opt) & (E_sorted < Eb2_opt)
m3 = (E_sorted >= Eb2_opt)

# # segment 1
# E1 = np.logspace(np.log10(E.min()), np.log10(Eb1_opt), 100)
# ccdf1 = (E1/E1[0])**(-(a1-1))

# # segment 2
E2 = np.logspace(np.log10(E.min()), np.log10(Eb2_opt), 100)
ccdf2 = (E2/E2[0])**(-(a2-1))

# # segment 3
E3 = np.logspace(np.log10(Eb2_opt), np.log10(E.max()), 100)
ccdf3 = ccdf2[-1]*(E3/Eb2_opt)**(-(a3-1))


# plt.loglog(E1,ccdf1,'-', label=f'α1={a1:.2f}')
plt.loglog(E2,ccdf2,'-', label=f'α2={a2:.2f}')
plt.loglog(E3,ccdf3,'-', label=f'α3={a3:.2f}')
plt.axvline(Eb1_opt, color='k', linestyle='--')
plt.axvline(Eb2_opt, color='k', linestyle='--')

plt.xlabel('E')
plt.ylabel('CCDF')
plt.legend()
plt.grid(True, which='both', alpha=0.3)
plt.savefig('pow_fit_3_test.png',dpi=300)
plt.show()