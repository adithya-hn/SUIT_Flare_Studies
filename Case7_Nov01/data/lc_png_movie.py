import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import glob
import matplotlib.animation as animation
import matplotlib.image as mpimg
import numpy as np
from datetime import datetime
import matplotlib.dates as mdates

path_dir='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case7_Nov01/data/aligned_crop_pngs/'
# Replace 'my_image.png' with the actual path to your image file

imgs=sorted(glob.glob(path_dir+'*.jpg'))

data1=(np.loadtxt(f'/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case7_Nov01/light_curve/csv_files/c7_NB04_lc_data.csv',skiprows=1,delimiter=',',dtype='str')).transpose() 


m_cls=datetime.fromisoformat('2024-11-13T14:18:00.000')
m_cls_p=datetime.fromisoformat('2024-11-13T14:31:00.000')

date_array1=data1[0]
# date_array2=data2[0]
time_array1=[]
time_array2=[]

for i in range(len(date_array1)):
    parsed_time = datetime.fromisoformat(date_array1[i])
    time_array1.append(parsed_time)

# for i in range(len(date_array2)):
#     parsed_time = datetime.fromisoformat(date_array2[i])
#     time_array2.append(parsed_time)

date=time_array1[0].strftime('%Y-%m-%d')
float_array1 = [float(string)  for string in data1[1]]
# float_array2 = [float(string)  for string in data2[1]]

# --- Set up plot ---
#fig, (ax_lc, ax_img) = plt.subplots(1, 2, figsize=(15, 6))
fig, (ax_img, ax_lc) = plt.subplots(2, 1, figsize=(8, 12), gridspec_kw={'height_ratios': [1.7, 1.2]})

# ax2=ax_lc.twinx()
#ax_lc.set_aspect(1.5) 
# Plot light curve
line_lc, = ax_lc.plot(time_array1, float_array1,'ko-',markersize=3, label="Mg II h")
# line_lc2, = ax2.plot(time_array2, float_array2,'bo-', markersize=3,label="Ca II h")
vline = ax_lc.axvline(time_array1[0], color='red', linestyle='--', label='Current Time')
ax_lc.set_ylabel("Mg II h counts",fontsize=12)
ax_lc.set_xlabel("Time (UT)",fontsize=12)
# ax_lc.set_yscale('log')
ax_lc.set_title("Light Curve")

# Load first image
img = mpimg.imread(imgs[0])
img_disp = ax_img.imshow(img)
ax_img.axis('off')
# ax_img.set_title(f"Image at {time_array1[0]}")

# --- Update function ---
def update(i):
    vline.set_xdata([time_array1[i]])

    img = mpimg.imread(imgs[i])
    img_disp.set_data(img)
    # ax_img.set_title(f"Image at {time_array1[i].strftime('%Y-%m-%d %H:%M:%S')}")

    return vline, img_disp

# --- Animate ---
ani = animation.FuncAnimation(fig, update, frames=len(time_array1), interval=300, blit=False)
time_formatter = mdates.DateFormatter('%H:%M') 
plt.gca().xaxis.set_major_formatter(time_formatter)
plt.figlegend(bbox_to_anchor=(0.001, 0.0015, 0.35, 0.5))
plt.tight_layout()
ani.save("lc_suit_aia_image.gif", writer='ffmpeg',dpi=260, fps=5)
plt.show()