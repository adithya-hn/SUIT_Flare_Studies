import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import glob
import matplotlib.animation as animation
import matplotlib.image as mpimg
import numpy as np
from datetime import datetime
import matplotlib.dates as mdates

path_dir='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case8_Nov01/localize_brightenings/Enhancment/NB04/'
# Replace 'my_image.png' with the actual path to your image file

imgs=sorted(glob.glob(path_dir+'*.png'))

data1=(np.loadtxt(f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case8_Nov01/light_curve/csv_files/c8_NB04_lc_data.csv',delimiter=',',skiprows=1,dtype='str')).transpose() 
#data2=(np.loadtxt(f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case8_Nov01/OVER_PLOT_CONTOURS/NB08_c8_lc_data.csv',delimiter=',',dtype='str')).transpose() 


m_cls=datetime.fromisoformat('2024-07-10T14:18:00.000')
m_cls_p=datetime.fromisoformat('2024-07-10T14:31:00.000')

date_array1=data1[0]
#date_array2=data2[0]
time_array1=np.array(data1[0], dtype='datetime64[us]')
#time_array2=[]

'''for i in range(len(date_array1)):
    parsed_time = datetime.fromisoformat(date_array1[i])
    time_array1.append(parsed_time)'''

# for i in range(len(date_array2)):
#     parsed_time = datetime.fromisoformat(date_array2[i])
#     time_array2.append(parsed_time)

date=time_array1[0]
float_array1 = [float(string)  for string in data1[1]]
#float_array2 = [float(string)  for string in data2[1]]

# --- Set up plot ---
#fig, (ax_lc, ax_img) = plt.subplots(1, 2, figsize=(15, 6))
fig, (ax_lc, ax_img) = plt.subplots(1, 2, figsize=(15, 6), gridspec_kw={'width_ratios': [1.8, 1.2]})

ax2=ax_lc.twinx()
#ax_lc.set_aspect(1.5) 
# Plot light curve
line_lc, = ax_lc.plot(time_array1, float_array1,'ko-',markersize=1, label="Mg II h")
#line_lc2, = ax2.plot(time_array2, float_array2,'bo-', markersize=1,label="Ca II h")
vline = ax_lc.axvline(time_array1[0], color='red', linestyle='--', label='Current Time')
ax_lc.set_ylabel("Counts")
#ax_lc.legend()
ax_lc.set_title("Light Curve")

# Load first image
img = mpimg.imread(imgs[0])
img_disp = ax_img.imshow(img)
ax_img.axis('off')
ax_img.set_title(f"SUIT Mg II h ")

# --- Update function ---
def update(i):
    vline.set_xdata([time_array1[i]])

    img = mpimg.imread(imgs[i])
    img = img[:-80, 80:]
    img_disp.set_data(img)
    #ax_img.set_title(f"Image at {time_array1[i]}")

    return vline, img_disp

# --- Animate ---
ani = animation.FuncAnimation(fig, update, frames=len(time_array1), interval=300, blit=False)
#plt.gca().xaxis.set_major_formatter(time_formatter)
time_formatter = mdates.DateFormatter('%H:%M')  # Format as HH:MM
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.figlegend(bbox_to_anchor=(0.001, 0.35, 0.35, 0.5))
plt.tight_layout()
ani.save("lightcurve_with_sliding_image.gif", writer='ffmpeg',dpi=300, fps=5)
plt.show()