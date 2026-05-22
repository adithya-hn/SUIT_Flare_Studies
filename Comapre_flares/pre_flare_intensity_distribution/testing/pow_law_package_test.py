"""
@Author      : Adithya H N
@Created On  : 2026-04-10
@Last Updated: 2026-04-10
@Project     : Pre-flare study
@Version     : 1.0

@Description
-----------
Brief description: checking the powerlaw fit from python package 'powerlaw'



"""


import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import powerlaw

data = np.loadtxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Comapre_flares/pre_flare_intensity_distribution/Flux_distrib/Peak_Flux.csv')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

E = np.sort(data)
N = len(E)
# ccdf_all = np.arange(N, 0, -1) / N

fit = powerlaw.Fit(E)
print(fit.power_law.alpha)
print(fit.power_law.xmin)

R, p = fit.distribution_compare('power_law', 'lognormal')
print(R,p)

import matplotlib.pyplot as plt
fig, ax = plt.subplots()

fit.plot_pdf(ax=ax, label='PDF')
fit.power_law.plot_pdf(ax=ax, label='Power law fit')

plt.legend()
plt.show()