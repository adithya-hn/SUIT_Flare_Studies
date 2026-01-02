
'''
Created on 14 oct 2025

'''


import glob
import sunpy.map
import astropy.units as u
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from sunkit_image.coalignment import mapsequence_coalign_by_match_template as mc_coalign

#-------------in put path------------------------

files = sorted(glob.glob("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case8_Nov01/data/1600_aligned/*NB06.fits"))

#-------------------------------------------------
os.makedirs("aligned_crop_", exist_ok=True)
os.makedirs("aligned_crop_pngs_", exist_ok=True)
os.makedirs("aligned_1600_pngs", exist_ok=True)
print(f'Number of files: {len(files)}')
'''

trX=[]
trY=[]
blX=[]
blY=[]
dates=[]

for i in range(len(files)):
    img = sunpy.map.Map(files[i])
    data = img.data
    nonzero = np.argwhere(data > 100)  # Find nonzero pixels
    ymin, xmin = nonzero.min(axis=0) # Get bounding box
    ymax, xmax = nonzero.max(axis=0)
    trX.append(xmax)
    trY.append(ymax)
    blX.append(xmin)
    blY.append(ymin)
    dates.append(img.date.datetime)
    # plt.imshow(data,origin='lower')
    # plt.scatter([xmin], [ymin], c='red', label='Top-right')
    # plt.scatter([xmax], [ymax], c='blue', label='Bottom-right')
    # plt.legend()
    # plt.savefig(f'aligned_1600_pngs/{img.meta["F_NAME"][:-4]}.png',dpi=200)
    # plt.close()

fig,ax1=plt.subplots(1,1, figsize=(10,5))
ax2=ax1.twinx()


ax1.plot(dates, trX,color='C0', label='Top Right X')
ax2.plot(dates, trY,color='C1',  label='Top Right Y')
ax1.plot(dates, blX,color='C2',  label='Bottom Left X')
ax2.plot(dates, blY,color='C3',  label='Bottom Left Y')

plt.xlabel('Time')
plt.ylabel('Pixel Coordinates')
plt.title('Stability of Image Alignment Over Time')
plt.figlegend(bbox_to_anchor=(0.001, 0.38, 0.4, 0.5))
plt.xticks(rotation=45)
plt.savefig('alignment_stability.png', dpi=200)
plt.close()

print(f'Crop coordinates: x1 {max(blX)}, y1 {max(blY)}, x2 { min(trX)}, y2 { min(trY)}')
'''
# x1,y1,x2,y2=85,125,620,760
#  #modified for cont
x1,y1,x2,y2=250,130,730,600

#os._exit(0) #------------------
seq=[]
for f in files:
    m = sunpy.map.Map(f)
    cropped = m.submap(bottom_left = [x1,y1]*u.pix,top_right = [x2,y2]*u.pix)
    seq.append(cropped)

Sequence_ = sunpy.map.MapSequence(seq)
aligned_maps=mc_coalign(Sequence_,layer_index=0,clip=True)
aligned_maps.peek()
ref_img_=aligned_maps[0]
check=[]
for j in range(len(aligned_maps)):  #account for incerted image, j=0 is refference image
    aligned_img=sunpy.map.Map(aligned_maps[j].data,aligned_maps[j].fits_header)
    aligned_img.meta['CRPIX1']=(ref_img_.meta['CRPIX1'])
    aligned_img.meta['CRPIX2']=(ref_img_.meta['CRPIX2'])
    aligned_img.meta['CRVAL1']=(ref_img_.meta['CRVAL1'])
    aligned_img.meta['CRVAL2']=(ref_img_.meta['CRVAL2'])
    fname= aligned_maps[j].meta.get('F_NAME')
    aligned_img.save('aligned_crop_/'+fname,overwrite=True) #need this for alignement refference
    check.append(aligned_img)
    
    fl_nm='aligned_crop_pngs_/'+fname[:-4]+'jpg'
    #Title_=aligned_img.meta.get('FTR_NAME') +' Filter: '+aligned_img.date
    fig=plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection=aligned_img)
    aligned_img.plot(title=False,axes=ax)
    plt.figtext(0.04, 0.04,aligned_img.date, color='white',weight='bold', fontsize=18)
    plt.draw()
    #plt.axis('off')
    plt.tight_layout()
    plt.colorbar()
    plt.savefig(fl_nm)
    plt.close()



    # fig=plt.figure()
    # ax=fig.add_subplot(111,projection=cropped)
    # cropped.plot()
    # plt.savefig(f'aligned_crop_pngs/{cropped.meta["F_NAME"][:-4]}.png',dpi=200)
    # plt.close()
    # cropped.save(f"aligned_crop/{cropped.meta["F_NAME"]}",overwrite=True)

