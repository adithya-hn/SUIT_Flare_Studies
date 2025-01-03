
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import scienceplots
import seaborn as sns
sns.set()
sns.set_style('ticks')
sns.set_context('notebook',font_scale=0.3)
plt.style.use('science')
plt.style.use(['science'])

plt.figure()

data=(np.loadtxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/June01_Flare/Difference_Image/Diff_img_data_NB08.dat',dtype='str')).transpose()

print(data.shape)
date_array=data[0]
time_array=[]
for i in range(len(data[0])):
    parsed_time = datetime.fromisoformat(date_array[i])
    time_array.append(parsed_time)
box_array = [float(string) for string in data[1]]
full_img_array = [float(string) for string in data[2]]
err_array = [float(string) for string in data[3]]

#plt.errorbar(time_array,box_array,yerr=err_array,fmt='ko--',capsize=0.5,markersize=0.5,linewidth=0.5)
plt.plot(time_array,box_array,'ko--',markersize=0.5,linewidth=0.5)
plt.xlabel('Time')
plt.ylabel('Counts')
#plt.legend(loc='best',fontsize=4,frameon=True)
plt.title('Difference image of box')
plt.savefig(f'NB08_Difference_box_image.png',dpi=800)
plt.show()