
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
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from astropy.wcs import WCS



flt='NB03'
nb3_maps=glob.glob(f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case3_June02/data/1600_aligned/*{flt}.fits')
m_seq=sunpy.map.Map(nb3_maps,sequence=True)


count=[]
date=[]
hist_data=[]
pathlib.Path(f"Diff_imgs/{flt}").mkdir(parents=True, exist_ok=True)
pathlib.Path("Diff_hists").mkdir(parents=True, exist_ok=True)

p1,p2,p3,p4=40,600,30,620
base_img_=np.zeros_like(m_seq[0].data[p1:p2,p3:p4])

for j in range(5):
    cutout=m_seq[j].data[p1:p2,p3:p4]
    #print((cutout.data).shape)
    base_img_+=(cutout*1000/m_seq[j].meta.get('CMD_EXPT'))

base_img=base_img_/5
base_map=sunpy.map.Map(base_img,m_seq[0].meta)
sigm_set=[]
mean_set=[]
for k in range(5):

    diff_img=((m_seq[k]).data[p1:p2,p3:p4]*1000/m_seq[k].meta.get('CMD_EXPT'))-(base_map.data)
    hist_data=diff_img.flatten()
    mean_val=np.mean(hist_data)
    sigma=np.std(hist_data)
    #print(sigma)
    sigm_set.append(sigma)
    mean_set.append(mean_val)

Thresh=np.mean(mean_set)+3*np.mean(sigm_set)
#print('Mean sigma: ',np.mean(sigm_set))

def gaussian(x,mean,amp,sigma):
    return  amp * np.exp(-(x - mean)**2 / (2 * sigma**2))
    

for i in range(len(m_seq)-1):
    fig = plt.figure()
    cutout=m_seq[i].data[p1:p2,p3:p4]
    diff_img=(cutout*1000/m_seq[i].meta.get('CMD_EXPT'))-(base_map.data)
    hist_data=diff_img.flatten()
    diff_img=diff_img
    mean_val=np.mean(hist_data)
    sigma=np.std(hist_data)
    Thresh=mean_val+3*sigma
    print(Thresh)
    
    Bins=np.arange(-1000,1000,10)
    
    img=np.where(diff_img>Thresh,diff_img,0)
    #img[:,300:350]=0
    #IMG=m_seq[i].data*1000/m_seq[i].meta.get('CMD_EXPT')
    
    plt.imshow(img,origin='lower',vmin=Thresh,vmax=300)
    #plt.imshow(m_seq[i+1].data[50:650,50:600],origin='lower',vmin=100,vmax=300,alpha=0.2)
    #plt.imshow(img)
    plt.colorbar(extend='both')
    plt.savefig(f'Diff_imgs/{flt}/{m_seq[i].date}.png')
    plt.close()
    count.append(np.sum(img))
    date.append(m_seq[i].date.datetime)

    '''
    n, bins, patches = plt.hist(hist_data, bins=Bins, color='gray', alpha=0.7, edgecolor='black')
    bin_centers = (bins[:-1] + bins[1:]) / 2
    params,cov=curve_fit(gaussian, bin_centers, n,p0=[0,6000,500])
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
#plt.yscale('log')
np.savetxt(f'{flt}_brightenings.csv',np.c_[date,count],delimiter=',',comments='',header='DateTime,Diff_count',fmt='%s')
plt.savefig(f'{flt}_base_simga_thresh_difference.png',dpi=300)
plt.show()

