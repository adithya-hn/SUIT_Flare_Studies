"""
@Author      : Adithya H N
@Created On  : 2026-04-13
@Last Updated: 2026-04-13
@Project     : pre-flare transients
@Version     : 1.0

@Description
-----------
Brief description: synthetic power law variation study
"""
import numpy as np
from scipy.optimize import minimize


def sample_truncated_powerlaw(alpha, xmin, xmax, size):
    """
    Inverse transform sampling for truncated power law.
    """
    u = np.random.uniform(size=size)

    a = xmin**(1-alpha)
    b = xmax**(1-alpha)

    return (a + u*(b-a))**(1/(1-alpha))


def sample_broken_powerlaw(alpha1, alpha2, xb, xmin, xmax, N):

    # relative probability in each segment
    I1 = (xb**(1-alpha1) - xmin**(1-alpha1)) / (1-alpha1)

    I2 = xb**(alpha2-alpha1) * \
         (xmax**(1-alpha2) - xb**(1-alpha2)) / (1-alpha2)

    frac1 = I1 / (I1 + I2)

    N1 = int(N * frac1)
    N2 = N - N1

    part1 = sample_truncated_powerlaw(alpha1, xmin, xb, N1)
    part2 = sample_truncated_powerlaw(alpha2, xb, xmax, N2)

    return np.sort(np.concatenate([part1, part2]))


def neg_logL_broken(params, data, xmin, xmax):

    alpha1, alpha2, xb = params

    if alpha1 <= 1 or alpha2 <= 1:
        return np.inf

    if not (xmin < xb < xmax):
        return np.inf

    data1 = data[data < xb]
    data2 = data[data >= xb]

    if len(data1) < 5 or len(data2) < 5:
        return np.inf

    N = len(data)

    norm = ((alpha1 - 1) * (alpha2 - 1)) / (
        xb * (
            alpha1 - alpha2
            + (alpha2 - 1) * (xmin / xb) ** (1 - alpha1)
            - (alpha1 - 1) * (xmax / xb) ** (1 - alpha2)
        )
    )

    if norm <= 0:
        return np.inf

    term1 = alpha1 * np.sum(np.log(data1 / xb))
    term2 = alpha2 * np.sum(np.log(data2 / xb))

    return -(N * np.log(norm) - term1 - term2)


def fit_broken_powerlaw(data, xmin, xmax):
    initial_guess = [1.7, 2.5, 1]

    bounds = [(0, 3), (1.01, 5), (0.8,1.1)]

    res = minimize(
        neg_logL_broken,
        initial_guess,
        args=(data, xmin, xmax),
        bounds=bounds,
        method='Nelder-Mead' )
    print(res.success, res.message)
    a1,a2,xb=res.x
    print(len(data[data<xb]), len(data[data>=xb]))
    return res.x  # alpha1, alpha2, xb
#-------------------------------------------------

true_alpha1 = 1.8
true_alpha2 = 2.6
true_xb     = 1.0

xmin = 0.36
xmax = 2.4

data_test = sample_broken_powerlaw(
    true_alpha1,
    true_alpha2,
    true_xb,
    xmin,
    xmax,
    N=5000
)

fit = fit_broken_powerlaw(data_test, xmin, xmax)

print("TRUE:")
print(true_alpha1, true_alpha2, true_xb)

print("\nRECOVERED:")
print(fit)

for i in range(10):
    d = sample_broken_powerlaw(1.8, 3.0, 1.0, xmin, xmax, 3000)
    fit = fit_broken_powerlaw(d, xmin, xmax)
    print(fit)

alpha1,alpha2,xb=fit

import matplotlib.pyplot as plt

bins = np.logspace(np.log10(xmin), np.log10(xmax), 30)

plt.hist(data_test, bins=bins, density=True, histtype='step')
plt.xscale('log')
plt.yscale('log')
plt.xlabel("x")
plt.ylabel("PDF")
plt.show()

E=data_test
def empirical_ccdf(E):
    n = len(E)
    return np.arange(n, 0, -1) / n

def model_ccdf(E, Emin, alpha):
    return (E / Emin) ** (-(alpha - 1))

eb2=xb
seg1 = E[E < eb2]
seg2 = E[E >= eb2]

model_ccdf_seg1 =  (seg1 / xb) ** (-(alpha1 - 1))
model_ccdf_seg2 =  (seg2 / xb) ** (-(alpha2 - 1))

plt.figure(figsize=(12,8))
plt.title('CCDF of all pre-flare transients',fontsize=24)

# Power-law fit (only above Emin)
plt.loglog(E, empirical_ccdf(E), '+',color='k', label='SUIT transeint events')
plt.loglog(seg1, model_ccdf_seg1, 'o',linewidth=2,label=fr'Power-law fit (alpha = {alpha1:.2f} )')
plt.loglog(seg2, model_ccdf_seg2, 'o',linewidth=2,label=fr'Power-law fit (alpha = {alpha2:.2f} )')
plt.show()


