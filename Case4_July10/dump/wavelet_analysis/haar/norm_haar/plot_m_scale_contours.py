
import numpy as np
import pywt
import glob
from astropy.io import fits
from tqdm import tqdm
import scipy.ndimage as ndi
from astropy.time import Time
import matplotlib.pyplot as plt
from skimage import filters, measure
from skimage.measure import label, regionprops
from skimage.morphology import disk, closing, opening,remove_small_objects
import seaborn as sns
import sunpy.map
import os
import pathlib
sns_cl3=sns.color_palette('bright')

fol_nm=os.getcwd() #Custom folder to save contour images
os.makedirs(fol_nm+'/comb_conts', exist_ok=True)
os.makedirs(fol_nm+'/m_scal_conts', exist_ok=True)

data  = fits.getdata('NB04_cube.fits')
data1 = fits.getdata('NB04_1_norm_haar_cube.fits')
data2 = fits.getdata('NB04_2_norm_haar_cube.fits')
data3 = fits.getdata('NB04_3_norm_haar_cube.fits')
data4 = fits.getdata('NB04_4_norm_haar_cube.fits')
data5 = fits.getdata('NB04_5_norm_haar_cube.fits')
data6 = fits.getdata('NB04_6_norm_haar_cube.fits')
data7 = fits.getdata('NB04_7_norm_haar_cube.fits')
data8 = fits.getdata('NB04_8_norm_haar_cube.fits')
nt,ny,nx=data1.shape

for i in range(nt):
    fig=plt.figure()
    plt.imshow(data[i,:,:], origin='lower', cmap='suit_nb04')
    comb   = data1[i,:,:]+data2[i,:,:]+data3[i,:,:]+data4[i,:,:]+data5[i,:,:]+data6[i,:,:]
    label1 = measure.label(data1[i,:,:],connectivity=2)
    label2 = measure.label(data2[i,:,:],connectivity=2)
    label3 = measure.label(data3[i,:,:],connectivity=2)
    label4 = measure.label(data4[i,:,:],connectivity=2)
    label5 = measure.label(data5[i,:,:],connectivity=2)
    label6 = measure.label(data6[i,:,:],connectivity=2)
    label7 = measure.label(data7[i,:,:],connectivity=2)
    label8 = measure.label(data8[i,:,:],connectivity=2)
    c_label= measure.label(comb,connectivity=2)
    if label1.max()>1:
        cleaned1 = remove_small_objects(label1, min_size=16)
    else:
        cleaned1=label1
    
    if label2.max()>1:
        cleaned2 = remove_small_objects(label2, min_size=16)
    else:
        cleaned2=label2

    if label3.max()>1:
        cleaned3 = remove_small_objects(label3, min_size=16)
    else:
        cleaned3=label3

    if label4.max()>1:
        cleaned4 = remove_small_objects(label4, min_size=16)
    else:
        cleaned4=label4

    if label5.max()>1:
        cleaned5 = remove_small_objects(label5, min_size=16)
    else:
        cleaned5=label5

    if label6.max()>1:
        cleaned6 = remove_small_objects(label6, min_size=16)
    else:
        cleaned6=label6
    
    if label7.max()>1:
        cleaned7 = remove_small_objects(label6, min_size=16)
    else:
        cleaned7=label7

    if label8.max()>1:
        cleaned8 = remove_small_objects(label6, min_size=16)
    else:
        cleaned8=label8

    if c_label.max()>1:
        cleaned = remove_small_objects(c_label, min_size=16)
    else:
        cleaned=c_label

    cont1 = measure.find_contours(cleaned1 > 0, level=0.5)
    cont2 = measure.find_contours(cleaned2 > 0, level=0.5)
    cont3 = measure.find_contours(cleaned3 > 0, level=0.5)
    cont4 = measure.find_contours(cleaned4 > 0, level=0.5)
    cont5 = measure.find_contours(cleaned5 > 0, level=0.5)
    cont6 = measure.find_contours(cleaned6 > 0, level=0.5)
    cont7 = measure.find_contours(cleaned7 > 0, level=0.5)
    cont8 = measure.find_contours(cleaned8 > 0, level=0.5)
    conts = measure.find_contours(cleaned > 0, level=0.5)
    
    for cont in cont1:
        plt.plot(cont[:, 1], cont[:, 0], linewidth=0.4, color=sns_cl3[1])
    
    for cont in cont2:
        plt.plot(cont[:, 1], cont[:, 0], linewidth=0.4, color=sns_cl3[2])
    
    for cont in cont3:
        plt.plot(cont[:, 1], cont[:, 0], linewidth=0.4, color=sns_cl3[3])
    
    for cont in cont4:
        plt.plot(cont[:, 1], cont[:, 0], linewidth=0.4, color=sns_cl3[4])
    
    for cont in cont5:
        plt.plot(cont[:, 1], cont[:, 0], linewidth=0.4, color=sns_cl3[5])
    
    for cont in cont6:
        plt.plot(cont[:, 1], cont[:, 0], linewidth=0.4, color=sns_cl3[6])

    for cont in cont7:
        plt.plot(cont[:, 1], cont[:, 0], linewidth=0.4, color=sns_cl3[6])
    
    for cont in cont8:
        plt.plot(cont[:, 1], cont[:, 0], linewidth=0.4, color=sns_cl3[6])

    plt.title(f'frame {i}')
    plt.savefig(f'm_scal_conts/img{i:02d}.png',dpi=200)
    plt.close()

    plt.figure()
    plt.title(f'frame{i}')
    plt.imshow(data[i,:,:], origin='lower', cmap='suit_nb04')
    for cont in conts:
        plt.plot(cont[:, 1], cont[:, 0], linewidth=0.4, color=sns_cl3[9])
    plt.savefig(f'comb_conts/img{i:02d}.png',dpi=200)
    plt.close()
