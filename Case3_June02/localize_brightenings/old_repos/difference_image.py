
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
'''
ax = fig.add_subplot(projection=m_seq.maps[0])
ani = m_seq.plot(axes=ax, norm=ImageNormalize(vmin=0, vmax=40000, stretch=SqrtStretch()))

plt.show()

m_seq_base = sunpy.map.Map([m - m_seq[0].quantity for m in m_seq[1:]], sequence=True)
m_seq_running = sunpy.map.Map(
    [m - prev_m.quantity for m, prev_m in zip(m_seq[1:], m_seq[:-1])],
    sequence=True
)

fig = plt.figure()
ax = fig.add_subplot(projection=m_seq_base.maps[0])
ani = m_seq_base.plot(axes=ax, title='Base Difference', norm=colors.Normalize(vmin=100, vmax=200), cmap='Greys_r')
plt.colorbar(extend='both', label=m_seq_base[0].unit.to_string())
plt.show()


fig = plt.figure()
ax = fig.add_subplot(projection=m_seq_running.maps[0])
ani = m_seq_running.plot(axes=ax, title='Running Difference', norm=colors.Normalize(vmin=100, vmax=200), cmap='Greys_r')'''
for i in range(len(m_seq)-1):
    fig = plt.figure()
    diff_img=m_seq[i+1].data-m_seq[i].data
    diff_img=diff_img[50:650,50:600]
    img=np.where(diff_img>3000,diff_img,0)
    
    plt.imshow(img,origin='lower')#,vmin=1000,vmax=2000)
    #plt.imshow(m_seq[i+1].data[50:650,50:600],origin='lower',vmin=1000,vmax=30000,alpha=0.1)
    #plt.imshow(img)
    plt.colorbar(extend='both')
    plt.savefig(f'Diff_imgs/{m_seq[i].date}.png')
    plt.close()
    count.append(np.sum(img))
    date.append(m_seq[i].date.datetime)

fig = plt.figure(figsize=(12,8))
plt.plot(date,count)
plt.savefig('running_th_difference.png',dpi=300)
plt.show()