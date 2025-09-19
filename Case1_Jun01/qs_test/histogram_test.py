import os
import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
import glob
import datetime

import pathlib
from astropy.coordinates import SkyCoord
import numpy as np
from scipy.ndimage import shift
from scipy.ndimage import gaussian_filter
from scipy.optimize import curve_fit

#Threshold values:

nb3T=11000
nb4T=11500
nb8T=3900
nb6T=86000
nb7T=295000
nb3Mx=14000
nb4Mx=15000
nb8Mx=4300
nb6Mx=95000
nb7Mx=310000
Filters=['NB03']
qs=[]
fit_qs=[]
time=[]
for fltr in Filters:
    plot_data=[]
    tot_count=[]
    dates=[]

    search_fold=f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case1_Jun01/data/cropped/crop_fits/{fltr}/' #Custom Folder

    
    print(f'Searching for {fltr} images in {search_fold} folder')
    fdir =search_fold 
    files = glob.glob(fdir + '*3'+fltr+'.fits')
    files=sorted(files, key=lambda file_name: datetime.datetime.strptime(os.path.basename(file_name).split('_')[5], "%Y-%m-%dT%H.%M.%S.%f"))
    files=files
    print('Total files:',len(files))

    ref_img=sunpy.map.Map(files[0])
    #base_fold='/home1/Data/Adithya/POC_Works/Jitter/'
    fl_date = datetime.datetime.fromisoformat(str(ref_img.date))
    fol_nm=os.getcwd() #str(fl_date.day).zfill(2)+'_'+str(fl_date.month).zfill(2)+'_'+str(fl_date.year).zfill(2)
    print(fol_nm)

    for i in range(len(files)):
        suitMap=sunpy.map.Map(files[i])
        img_head=suitMap.fits_header
        suit_data=suitMap.data#gaussian_filter(suitMap.data,sigma=sigma)
        #print('before norm:',suit_data[100,100])
        alned_data=suit_data*1000/int(suitMap.meta.get('CMD_EXPT'))
        
        #print('after norm:',alned_data[100,100],'Expos:',int(suitMap.meta.get('MEAS_EXP')))

        suit_map=sunpy.map.Map(alned_data,img_head)

        Bins=np.arange(1000,20000,10)
        
        # Plot histogram
        plt.figure(figsize=(8,6))
        flt_data=suit_map.data.flatten()
        flt_data=flt_data[flt_data>1]
        n, bins, patches = plt.hist(flt_data, bins=Bins, color='gray', alpha=0.7, edgecolor='black')

        max_bin_index = np.argmax(n)
        qs_intensity = (bins[max_bin_index] + bins[max_bin_index+1]) / 2

        print(f"Estimated Quiet Sun (QS) intensity: {qs_intensity}")
        plt.legend([f'QS: {qs_intensity}'])
        #plt.axvline(qs_intensity, color='red', linestyle='dashed', linewidth=1)
        #plt.xscale('log')
        
        bin_centers = (bins[:-1] + bins[1:]) / 2
        qs_guess = bin_centers[max_bin_index]

        md = np.median(flt_data)

        sig_best=int(np.median(abs(flt_data-md))/0.6745)
        def gaussian(x, amp, mean, sigma):
            return amp * np.exp(-(x - mean)**2 / (2 * sigma**2))
        
        def triple_gaussian(x, amp1, mean1, sigma1, amp2, mean2, sigma2, amp3, mean3, sigma3):
            gauss1 = amp1 * np.exp(-(x - mean1)**2 / (2 * sigma1**2))
            gauss2 = amp2 * np.exp(-(x - mean2)**2 / (2 * sigma2**2))
            gauss3 = amp3 * np.exp(-(x - mean3)**2 / (2 * sigma3**2))
            return gauss1 + gauss2 + gauss3
        
        # Peak of histogram as QS guess
        max_bin_index = np.argmax(n)
        qs_guess = bin_centers[max_bin_index]

        # Rough guesses
        p0 = [
            np.max(n), qs_guess, 100,  # First Gaussian (QS)
            np.max(n)/5, qs_guess+500, 200,  # Second Gaussian (active region)
            np.max(n)/10, qs_guess+1000, 300 ] # Third Gaussian (brighter features)
        try:
            params, cov = curve_fit(triple_gaussian, bin_centers, n, p0=p0,maxfev=5000)

            amp1, mean1, sigma1, amp2, mean2, sigma2, amp3, mean3, sigma3 = params


            print(f"QS (Gaussian 1): Mean = {mean1:.2f}, Sigma = {sigma1:.2f}")
            print(f"Active Region (Gaussian 2): Mean = {mean2:.2f}, Sigma = {sigma2:.2f}")
            print(f"Bright Features (Gaussian 3): Mean = {mean3:.2f}, Sigma = {sigma3:.2f}")

            
            # Total fit
            x_fit = np.linspace(bin_centers.min(), bin_centers.max(), 2000)
            y_fit = triple_gaussian(x_fit, *params)
            plt.plot(x_fit, y_fit, 'r-', lw=2, label='Triple Gaussian Fit')

            # Individual Gaussians
            y_gauss1 = amp1 * np.exp(-(x_fit - mean1)**2 / (2 * sigma1**2))
            y_gauss2 = amp2 * np.exp(-(x_fit - mean2)**2 / (2 * sigma2**2))
            y_gauss3 = amp3 * np.exp(-(x_fit - mean3)**2 / (2 * sigma3**2))

            plt.plot(x_fit, y_gauss1, 'b--', label=f'Gaussian 1 (QS): μ={mean1:.2f}')
            plt.plot(x_fit, y_gauss2, 'g--', label=f'Gaussian 2: μ={mean2:.2f}')
            plt.plot(x_fit, y_gauss3, 'm--', label=f'Gaussian 3: μ={mean3:.2f}')

            plt.xlabel('Pixel Intensity')
            plt.ylabel('Number of Pixels')
            plt.title('Histogram of SUIT ROI with Triple Gaussian Fit')
            plt.legend()
            plt.grid(True)
            plt.axvline(qs_guess+0.5*sig_best, color='blue', linestyle=':')
            plt.savefig(f'histogram_test/Histogram_{fltr}_{i}.png', dpi=300, bbox_inches='tight')
            plt.close()# Estimate initial QS guess from histogram peak
            qs.append(qs_intensity)
            fit_qs.append(mean1)
            time.append(suit_map.date.datetime)

        except :
            print("Gaussian fit failed, skipping...")
            pass

np.savetxt(f'QS_values_{fltr}.csv',np.c_[time,qs],delimiter=',',fmt='%s',header='Date,QS_Intensity')
plt.figure(figsize=(10, 5))
plt.plot(time, qs, marker='o',markersize=0.5)
plt.plot(time, fit_qs, marker='o',markersize=0.5)
plt.gcf().autofmt_xdate()
plt.xlabel("Time")
plt.ylabel("Quiet Sun Intensity")
plt.title(f"Quiet Sun Intensity Variation {fltr} Å")
plt.grid(True)

plt.savefig(f'QS_values_{fltr}.png', dpi=300, bbox_inches='tight')
plt.show()
        # Find bin with maximum count (mode of histogram)
        


