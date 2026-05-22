

"""
@Author      : Adithya H N
@Created On  : 2026-04-12
@Last Updated: 2026-04-13
@Project     : Pre-flare study
@Version     : 1.0

@Description
-----------
Brief description: Applying MLE base powerlaw base on the paper 'Maximum likelihood estimation of the broken power law
spectral parameters with detector design applications by Leonard W. Howell
min threshold was set maually


"""


import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

data = np.loadtxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Comapre_flares/pre_flare_intensity_distribution/Flux_distrib/Peak_Flux.csv')
data=np.sort(data[data>0.33])


data = np.sort(data[data > 0.33])
EA   = 0.33           # lower flux bound (x_min)
EB   = data.max()     # upper flux bound (x_max)
N    = len(data)

def neg_logL_howell(params, data, EA, EB):

    a1, a2, Ek = params
    if a1 <= 1 or a2 <= 1:          return np.inf
    if a2 <= a1:                     return np.inf   # Howell condition
    if not (EA < Ek < EB):           return np.inf

    d1 = data[data < Ek]
    d2 = data[data >= Ek]
    if len(d1) < 20 or len(d2) < 20:  return 1e12

    # Normalization
    I1 = (Ek**(1 - a1) - EA**(1 - a1)) / (1 - a1)
    I2 = Ek**(a2 - a1) * (EB**(1 - a2) - Ek**(1 - a2)) / (1 - a2)
    inv_A = I1 + I2
    if inv_A <= 0 or not np.isfinite(inv_A): return 1e12

    # Log-likelihood 
    ll = (N * (-np.log(inv_A))
          - a1 * np.sum(np.log(d1))
          + (a2 - a1) * len(d2) * np.log(Ek)  
          - a2 * np.sum(np.log(d2)))

    return -ll if np.isfinite(ll) else 1e12


def fit_howell_broken_pl(data, xmin,xmax):
    best_l=np.inf
    a1_vals=np.arange(0,2,0.01)
    a2_vals=np.arange(2,5,0.01)
    eb_vals=np.arange(0.8,1.8,0.01)
    l_vals=[]

    for a1 in a1_vals:
        for a2 in a2_vals:
            for eb in eb_vals:
                param=[a1,a2,eb]
                neg_l=neg_logL_howell(param,data,xmin,xmax)
                l_vals.append(neg_l)
                if neg_l<best_l:
                    best_l=neg_l
                    best_param=[a1,a2,eb,best_l]
                    

    best_vals=min(l_vals)
    return best_param


# print(fit_howell_broken_pl(data, min(data),max(data)))``
alpha1, alpha2, Ek,max_L =fit_howell_broken_pl(data, min(data),max(data))
print(f"alpha1 = {alpha1:.4f}")
print(f"alpha2 = {alpha2:.4f}")
print(f"Ek     = {Ek:.4f}")
print(f"n1 = {np.sum(data < Ek)}, n2 = {np.sum(data >= Ek)}")
print(f"Maximum lilklihood: ",max_L)

#------------------------------------------
E=data
xb=Ek
# xb=1
# alpha1=0.1
# alpha2=2.8


def empirical_ccdf(E):
    n = len(E)
    return np.arange(n, 0, -1) / n

xgrid = np.logspace(np.log10(min(data)), np.log10(max(data)), 500)
# normalization denominator
I1 = (xb**(1-alpha1) - EA**(1-alpha1)) / (1-alpha1)
I2 = xb**(alpha2-alpha1) * (max(data)**(1-alpha2) - xb**(1-alpha2)) / (1-alpha2)

A = 1 / (I1 + I2)

ccdf_model = []

for x in xgrid:

    if x < xb:
        tail1 = (xb**(1-alpha1) - x**(1-alpha1)) / (1-alpha1)
        tail2 = xb**(alpha2-alpha1) * (max(data)**(1-alpha2) - xb**(1-alpha2)) / (1-alpha2)
        ccdf_model.append(A * (tail1 + tail2))

    else:
        tail2 = xb**(alpha2-alpha1) * (max(data)**(1-alpha2) - x**(1-alpha2)) / (1-alpha2)
        ccdf_model.append(A * tail2)

ccdf_model = np.array(ccdf_model)
plt.figure(figsize=(12,8))

plt.loglog(E, empirical_ccdf(E), '+', color='k', label='Data')
plt.loglog(xgrid, ccdf_model, '-', lw=3, label='Broken PL MLE')

plt.axvline(xb, ls='--', color='r', label=f'Break = {xb:.2f}')

plt.legend()
plt.show()


xmax=max(data)
xmin=min(data)
def neg_logL_single(alpha):
    if alpha <= 1: return 1e12
    norm = (1 - alpha) / (xmax**(1-alpha) - xmin**(1-alpha))
    return -(N * np.log(norm) - alpha * np.sum(np.log(data)))

r_single = minimize_scalar(neg_logL_single, bounds=(1.001, 5), method='bounded')

# LRT statistic (chi-sq with 2 dof)
lrt = 2 * (r_single.fun - neg_logL_howell([alpha1, alpha2, xb], data, xmin, xmax))
print(f"LRT stat = {lrt:.2f}  (p ~ chi2 cdf with 2 dof)")
delta_aic = 2*2 - lrt  # broken has 2 extra params
print(f"Delta AIC (broken - single) = {delta_aic:.2f}  (<0 favors broken)")

# Lower-branch slope guide
x1 = np.logspace(np.log10(EA), np.log10(xb), 100)
y1 = ccdf_model[np.argmin(np.abs(xgrid - EA))] * (x1 / EA)**(-(alpha1 - 1))

# Upper-branch slope guide
x2 = np.logspace(np.log10(xb), np.log10(EB), 100)

# Anchor upper slope at break for continuity
y_break = ccdf_model[np.argmin(np.abs(xgrid - xb))]
y2 = y_break * (x2 / xb)**(-(alpha2 - 1))

plt.loglog(x1, y1, '--', lw=2,
           label=fr'CCDF slope = $-(\alpha_1-1)={-(alpha1-1):.2f}$')

plt.loglog(x2, y2, '--', lw=2,
           label=fr'CCDF slope = $-(\alpha_2-1)={-(alpha2-1):.2f}$')
plt.show()

#----------------------------------
# plt.figure(figsize=(12,8))
# plt.title('CCDF of all pre-flare transients',fontsize=24)

# Power-law fit (only above Emin)
# plt.loglog(E, empirical_ccdf(E), '+',color='k', label='SUIT transeint events')
# plt.loglog(seg1, model_ccdf_seg1, 'o',linewidth=2,label=fr'Power-law fit (alpha = {alpha1:.2f} )')
# plt.loglog(seg2, model_ccdf_seg2, 'o',linewidth=2,label=fr'Power-law fit (alpha = {alpha2:.2f} )')

# plt.show()

# bins = np.logspace(np.log10(min(data)), np.log10(max(data)), 20)

# plt.hist(data, bins=bins, density=True, histtype='step')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel("x")
# plt.ylabel("PDF")
# plt.show()