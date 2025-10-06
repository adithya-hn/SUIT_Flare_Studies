
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import ImageNormalize, SqrtStretch
import tqdm
import sunpy.data.sample
import sunpy.map
import glob
import pathlib
import astropy.units as u
from pylab import *
import matplotlib.dates as mdates
from scipy.stats import mode
from scipy.optimize import curve_fit

nb3_maps=glob.glob('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/data/1600_aligned/*NB08.fits')
m_seq=sunpy.map.Map(nb3_maps,sequence=True)

count=[]
date=[]
hist_data=[]
pathlib.Path("Diff_imgs/NB08").mkdir(parents=True, exist_ok=True)
pathlib.Path("Diff_hists").mkdir(parents=True, exist_ok=True)

base_img=np.zeros_like(m_seq[0].data)
for j in range(5):
    base_img+=(m_seq[j].data*1000/m_seq[j].meta.get('CMD_EXPT'))
base_img=base_img/5
base_map=sunpy.map.Map(base_img,m_seq[0].meta)
diff_img=(m_seq[2].data*1000/m_seq[2].meta.get('CMD_EXPT'))-(base_map.data)
hist_data=diff_img.flatten()
mean_val=np.mean(hist_data)
sigma=np.std(hist_data)
print('Threshold: ',mean_val+3*sigma)

def gaussian(x,mean,amp,sigma):
    return  amp * np.exp(-(x - mean)**2 / (2 * sigma**2))
    

for i in range(len(m_seq)-1):
    fig = plt.figure()
    diff_img=(m_seq[i].data*1000/m_seq[i].meta.get('CMD_EXPT'))-(base_map.data)
    diff_img=diff_img[50:600,50:600]
    hist_data=diff_img.flatten()
    mean_val=np.mean(hist_data)
    sigma=np.std(hist_data)
    Thresh=mean_val+3*sigma
    Bins=np.arange(-1000,1000,10)
    '''
    n, bins, patches = plt.hist(hist_data, bins=Bins, color='gray', alpha=0.7, edgecolor='black')
    bin_centers = (bins[:-1] + bins[1:]) / 2
    params,cov=curve_fit(gaussian, bin_centers, n,p0=[0,6000,500])'''
    img=np.where(diff_img>Thresh,diff_img,0)
    #IMG=m_seq[i].data*1000/m_seq[i].meta.get('CMD_EXPT')
    plt.imshow(img,origin='lower',vmin=50,vmax=100)
    #plt.imshow(m_seq[i+1].data[50:650,50:600],origin='lower',vmin=100,vmax=300,alpha=0.2)
    #plt.imshow(img)
    plt.colorbar(extend='both')
    plt.savefig(f'Diff_imgs/NB08/{m_seq[i].date}.png')
    plt.close()
    count.append(np.sum(img))
    date.append(m_seq[i].date.datetime)
    
    '''
    # Total fit
    plt.bar(bin_centers, n, width=(bins[1]-bins[0]), color='gray', alpha=0.6, label='Histogram')
    x_fit = np.linspace(bin_centers.min(), bin_centers.max(), 2000)
    y_fit = gaussian(x_fit, *params)
    # Extract fitted values
    fit_vals = gaussian(bin_centers, *params)
    # Residuals
    residuals = n - fit_vals
    chi2 = np.sum((residuals**2) / (fit_vals + 1e-6))   # avoid divide by 0
    ndof = len(Bins) - len(params)  # degrees of freedom
    reduced_chi2 = chi2 / ndof
    print("Fit parameters:", params)
    print(f"Reduced chi-square: {reduced_chi2:.3f}")
    plt.plot(x_fit, y_fit, 'r-', lw=2, label='Gaussian Fit')
    #plt.axvline()
    plt.savefig(f'Diff_hists/Diff_Histogram{m_seq[i].date}.png',dpi=300)
    plt.close()'''
    


rc('axes', linewidth=1.2)
plt.rcParams["xtick.major.size"] = 10
fig,axs=plt.subplots(1,1, figsize=(10,5))
#ax2 = axs.twinx()
axs.xaxis.set_tick_params(size=0.5)
axs.yaxis.set_tick_params(size=0.5)
axs.tick_params(axis='both', direction='in', length=6, width=1)
axs.tick_params(which='minor', direction='in', length=3, width=1)
axs.yaxis.set_ticks_position('both')
axs.xaxis.set_ticks_position('both')
axs.minorticks_on()
axs.plot(date,count)
time_formatter = mdates.DateFormatter('%H:%M') 
plt.savefig('nb08base_simga_thresh_difference.png',dpi=300)
plt.show()

