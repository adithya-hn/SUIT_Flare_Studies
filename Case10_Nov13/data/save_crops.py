import glob
import sunpy.map
import astropy.units as u
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

os.makedirs("aligned_crop", exist_ok=True)
os.makedirs("aligned_crop_pngs", exist_ok=True)
files = sorted(glob.glob("/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case10_Nov13/data/1600_aligned/*.fits"))
x1, x2, y1, y2 = 110, 750, 110, 750

trX=[]
trY=[]
blX=[]
blY=[]
dates=[]

for i in range(len(files)):
    img = sunpy.map.Map(files[i])
    data = img.data
    nonzero = np.argwhere(data > 0)  # Find nonzero pixels
    ymin, xmin = nonzero.min(axis=0) # Get bounding box
    ymax, xmax = nonzero.max(axis=0)
    trX.append(xmax)
    trY.append(ymin)
    blX.append(xmin)
    blY.append(ymin)
    dates.append(img.date.datetime)

plt.plot(dates, trX, label='Top Right X')
plt.plot(dates, trY, label='Top Right Y')
plt.plot(dates, blX, label='Bottom Left X')
plt.plot(dates, blY, label='Bottom Left Y')
plt.xlabel('Time')
plt.ylabel('Pixel Coordinates')
plt.title('Stability of Image Alignment Over Time')
plt.legend()
plt.xticks(rotation=45)
plt.savefig('alignment_stability.png', dpi=200)
plt.close()

# Get bounding box
ymin, xmin = nonzero.min(axis=0)
ymax, xmax = nonzero.max(axis=0)

print(f"Top-right corner:    (x={xmax}, y={ymin})")

for f in files:
    m = sunpy.map.Map(f)
    cropped = m.submap(bottom_left = [x1,y1]*u.pix,top_right = [x2,y2]*u.pix)
    fig=plt.figure()
    ax=fig.add_subplot(111,projection=cropped)
    cropped.plot()
    #plt.imshow(cropped.data,origin='lower')
    plt.savefig(f'aligned_crop_pngs/{cropped.meta["F_NAME"][:-4]}.png',dpi=200)
    plt.close()


    cropped.save(f"aligned_crop/{cropped.meta["F_NAME"]}",overwrite=True)

