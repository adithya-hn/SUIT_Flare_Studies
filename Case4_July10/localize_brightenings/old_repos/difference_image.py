
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

nb3_maps=glob.glob('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case3_June02/data/1600_aligned/*NB03.fits')
m_seq=sunpy.map.Map(nb3_maps,sequence=True)

count=[]
date=[]
pathlib.Path("Diff_imgs").mkdir(parents=True, exist_ok=True)

for i in range(len(m_seq)-1):
    fig = plt.figure()
    diff_img=m_seq[i+1].data-m_seq[i].data
    diff_img=diff_img[50:650,50:600]
    img=np.where(diff_img>1000,diff_img,0)
        
    
    plt.imshow(img,origin='lower',vmin=1000,vmax=2000)
    plt.imshow(m_seq[i+1].data[50:650,50:600],origin='lower',vmin=1000,vmax=30000,alpha=0.1)
    #plt.imshow(img)
    plt.colorbar(extend='both')
    plt.savefig(f'Diff_imgs/{m_seq[i].date}.png')
    plt.close()
    count.append(np.sum(img))
    date.append(m_seq[i].date.datetime)

plt.plot(date,count)
plt.show()