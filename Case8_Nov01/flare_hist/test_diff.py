'''
Created on 01/10/2025
@author: adithya-hn

- To make differnce image of aligend suit images. 
'''




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
from skimage import filters, measure
from skimage.measure import label, regionprops
from skimage.morphology import disk, closing, opening,remove_small_objects

from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
import align_suit_fltr_to_aia


def diff_img(data_folder,flt,c1=0,c2=0,thresh_sig=4,clip_bar=False,get_histograms=False):
    flt_maps=glob.glob(data_folder)
    m_seq=sunpy.map.Map(flt_maps,sequence=True)


    count=[]
    count_er=[]
    date=[]
    hist_data=[]
    sigm_set=[]
    mean_set=[]
    img_count=[]

    pathlib.Path(f"Diff_imgs/{flt}").mkdir(parents=True, exist_ok=True)
    pathlib.Path(f"Enhancment/{flt}").mkdir(parents=True, exist_ok=True)

    '''
    base_img_=np.zeros_like(m_seq[0].data)
    base_img_er=np.zeros_like(m_seq[0].data)
    for j in range(5):
        cutout=m_seq[j].data
        #print((cutout.data).shape)
        base_img_+=(cutout*1000/m_seq[j].meta.get('CMD_EXPT'))
        base_img_er+=(np.sqrt(cutout)*1000/m_seq[j].meta.get('CMD_EXPT'))**2
    base_imgEr=np.sqrt(base_img_er)/5'''

    data_stack = np.stack([(m_seq[i].data*1000/m_seq[i].meta.get('CMD_EXPT'))for i in range(5)])
    base_img=np.median(data_stack, axis=0)

    data_stack_ = np.stack([(m_seq[i].data)for i in range(5)])
    base_img_er=np.sqrt(np.median(data_stack_, axis=0))*1000/m_seq[1].meta.get('CMD_EXPT')
    print([m_seq[i].meta.get('CMD_EXPT')for i in range(5)])

    # base_mad_sig=np.median(np.abs(data_stack-base_img),axis=0)/0.6745

    # plt.imshow(base_mad_sig,vmin=0,vmax=500)
    # plt.colorbar()
    # plt.savefig('Base_mad_sig.png',dpi=300)
    # plt.close()

    # plt.imshow(base_img_er_,vmin=0,vmax=500)
    # plt.colorbar()
    # plt.savefig('Base_poison_er.png',dpi=300)
    # plt.close()

    # plt.imshow(base_img_er_-base_mad_sig,cmap='hmimag',vmax=200,vmin=-200)
    # plt.colorbar()
    # plt.savefig('Basemapd_sig-poison.png',dpi=300)
    # plt.close()

    # print('Base er',np.sum(base_img_er),np.sum(base_mad_sig),np.sum(base_img_er-base_mad_sig))


    base_map=sunpy.map.Map(base_img,m_seq[0].meta) #header of any of 5, considering first one
    
    
    def double_gaussian(x, amp1, mean1, sigma1, amp2, mean2, sigma2):
            gauss1 = amp1 * np.exp(-(x - mean1)**2 / (2 * sigma1**2))
            gauss2 = amp2 * np.exp(-(x - mean2)**2 / (2 * sigma2**2))
            return gauss1 + gauss2

    def gaussian(x,mean,amp,sigma):
        return  amp * np.exp(-(x - mean)**2 / (2 * sigma**2))

    for i in range(len(m_seq)):
        fig = plt.figure()
        cutout=m_seq[i].data
        diff_img=(cutout*1000/m_seq[i].meta.get('CMD_EXPT'))-(base_map.data)
        diff_img_er=np.sqrt((np.sqrt(cutout)*1000/m_seq[i].meta.get('CMD_EXPT'))**2+(base_img_er)**2)

        hist_data=np.array(diff_img.flatten(),dtype='int')
        mean_val=np.mean(hist_data)
        median_val=np.median(hist_data)
        sigma=np.std(hist_data)
        sigm_set.append(sigma)
        mean_set.append(mean_val)
        mad_sig=np.median(np.abs(hist_data-np.median(hist_data)))/0.6745
        v_Thresh=median_val+thresh_sig*mad_sig
        #print('Mean-median',mean_val,median_val,'| Threshold valus: ',v_Thresh) 
        img=np.where(diff_img>v_Thresh,diff_img,0)

        #------------
        binary_image = diff_img > v_Thresh# True where pixel value > threshold
        labels = measure.label(binary_image,connectivity=2)
        if labels.max()>1:
           cleaned = remove_small_objects(labels, min_size=16)
        else:
           cleaned=labels
        #cleaned=labels
        mask=binary_image>0
        print('------------')
        #contours = measure.find_contours(binary_image, level=0.5)
        contours = measure.find_contours(cleaned > 0, level=0.5)
        fig=plt.figure()
        ax = fig.add_subplot(projection=m_seq[i])
        #fig,ax=plt.subplots()
        color_scale='gray'
        imax=25000
        if flt=='NB04':
            color_scale='suit_nb04'
            imax=28000
        if flt=='NB08':
            color_scale='suit_nb08'
            imax=10000
        if flt=='NB01':
            color_scale='suit_nb01'
            imax=11000
        if flt=='NB02':
            color_scale='suit_nb02'
            imax=35000
        if flt=='NB05':
            color_scale='suit_nb05'
            imax=48000
        if flt=='NB06':
            color_scale='suit_nb06'
            imax=95000
        if flt=='NB07':
            color_scale='suit_nb07'
            imax=270000
        if flt=='NB04':
            color_scale='suit_nb04'
       
        if flt=='NB03':
            color_scale='suit_nb03'
        # plt.imshow(m_seq[i].data, origin='lower', cmap=color_scale, vmin=1000, vmax=imax)
        # for contour in contours:
        #     ax.plot(contour[:, 1], contour[:, 0], linewidth=0.6, color='red')
        # ax.tick_params(labelsize=20) 
        # plt.colorbar(label='DN/s')
        # plt.title(f'{str(m_seq[i].date)[11:-4]}',fontsize=22)
        # plt.savefig(f'Enhancment/{flt}/{(m_seq[i].meta["F_NAME"])[:-4]}.png', dpi=150)
        # plt.close()
        

        th_diff_img=mask.astype(int)*diff_img
        th_diff_img_er= np.sqrt(np.sum((mask.astype(int)*diff_img_er)**2))
        th_img=mask.astype(int)*m_seq[i].data
        
        fig = plt.figure()
        plt.imshow(th_img)
        plt.savefig(f'{(m_seq[i].meta["F_NAME"])[:-4]}.png')
        plt.show()
        print('img-sum',np.sum(th_img))
        print('diff img-sum',np.sum(th_diff_img))

        


        # plt.imshow(img,origin='lower',vmin=Thresh,vmax=300)

        # plt.imshow(img,origin='lower',vmin=0,vmax=50)
        # plt.colorbar(extend='both')
        # plt.savefig(f'Diff_imgs/{flt}/{m_seq[i].date}.png')
        # plt.close()
        count.append(np.sum(th_diff_img))
        count_er.append(th_diff_img_er)
        img_count.append(np.sum(th_img))
        date.append(m_seq[i].date.datetime)
        
        
        
        
    print(len(date), len(mean_set), len(sigm_set), len(count), len(img_count))
    np.savetxt(f'Diff_img_data_{flt}.csv',np.c_[date,mean_set,sigm_set,count,count_er,img_count],delimiter=',',comments='',header='Date,Mean,Std,diff_count,Diff_error,Img_count',fmt='%s')


if __name__=='__main__':
    #--------------------Input params-------------------------------

    
    p1,p2,p3,p4=40,670,40,680
    flt='NB04' #NB04,NB08
    data_folder=f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/data/aligned_crop/*{flt}.fits'
    clip_bar=False
    c1=70
    c2=130
    thresh_sig=5
    get_histograms=False
    diff_img(data_folder,flt,c1,c2,thresh_sig,clip_bar)