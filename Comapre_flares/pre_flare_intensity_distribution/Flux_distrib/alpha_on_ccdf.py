import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

data = np.loadtxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Comapre_flares/pre_flare_intensity_distribution/Flux_distrib/Peak_Flux.csv')

data=np.sort(data)
E=data
def empirical_ccdf(E):
    n = len(E)
    return np.arange(n, 0, -1) / n

Eb1_opt=0.36
Eb2_opt=1.0
a2=1.71
a3=3.28
# segment 1
# E1 = np.logspace(np.log10(E.min()), np.log10(Eb1_opt), 100)
# ccdf1 = (E1/E1[0])**(-(a1-1))

# segment 2
E2 = np.logspace(np.log10(Eb1_opt), np.log10(Eb2_opt), 100)
ccdf2 = 0.65*(E2/Eb1_opt)**(-(a2-1)) #ccdf1[-1]*

# a3=3.1

# # segment 3
E3 = np.logspace(np.log10(Eb2_opt), np.log10(E.max()), 100)
ccdf3 = ccdf2[-1]*(E3/Eb2_opt)**(-(a3-1))
plt.loglog(data,empirical_ccdf(data), '+',color='k', label='SUIT transeint events')
plt.loglog(E2,ccdf2,'-',markersize=9 ,label=fr'$\alpha$ 2={a2:.2f}' ,alpha=0.8)
plt.loglog(E3,ccdf3,'-', label=fr'$\alpha$ 3={a3:.2f}',alpha=0.8)

plt.title('CCDF of pre-flare transients', fontsize=20)
plt.xlabel(r'Peak flux ($erg^-1 \times cm^-2$)', fontsize=18)
plt.ylabel('CCDF', fontsize=18)
plt.legend(fontsize=14)
plt.grid(True, which='both', alpha=0.3)
# plt.savefig('ccdf_mle_2bkn_pow_law.pdf',dpi=300)
plt.savefig('alpha_on_ccdf.png',dpi=300)
plt.show()