# File operations
import glob as glob
# Visualizations
import sunpy.map
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings('ignore')
# Include local library paths
import sys
sys.path.append('/home/adithya/Adithya_repos/pil_detection/utils')

# Import local libraries
from region_detection import pos_neg_detection
from pil_detection import detection
from video_loading import video_loader

# # Define HARP AR number as well as input and output path
# HARP_AR = 377
input_path = '/media/adithya/Adi_disk4/SUIT_flare_work/case5_jul10/data/HMI/HMI_cutouts/'
output_path = '/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case5_July10/mag/PIL_detection/PIL_images/'

# # Create class object by initializing path to samples and the HARP number of the sample input
dt = pos_neg_detection()
pil_dt = detection(input_path)
# loader = video_loader()

# # Check whether output path exists
pil_dt.check_outpath(output_path)


STRENGTH_FILTER = 100
BUFFER_SIZE = 4
PRESERVED_FLUX_RATIO = 0.95
MIN_MPIL_SIZE = 14

# Count the number of samples
data = glob.glob('/media/adithya/Adi_disk4/SUIT_flare_work/case5_jul10/data/HMI/HMI_cutouts/*.fits')
print('Total number of sample images: ' + str(len(data)))

# Specify the index of desired magnetogram patch
# Our sample includes 50 elements (indices 0 to 49)
data_index = 1

hmi_magmap = sunpy.map.Map(data[data_index])


pil_orig, label_orig = pil_dt.PIL_detect(pos_gauss = STRENGTH_FILTER, neg_gauss= -STRENGTH_FILTER, size_kernel = BUFFER_SIZE)
data_num = pil_dt.check_header(hmi_magmap)
masked_pil = dt.mask_pil(label_orig[data_num])
strength_label = pil_dt.filter_by_strength(threshold = PRESERVED_FLUX_RATIO)
masked_filter_RoPIs = dt.mask_pil(strength_label[data_num])
thin_df, thin_binary = pil_dt.thin_strength_label(strength_label)
masked_thinned_MPIL = dt.mask_pil(thin_binary[data_num])
pil_msk=sunpy.map.Map(masked_thinned_MPIL,hmi_magmap.meta)
        

fig = plt.figure(figsize=(8,6))
#hmi_magmap.plot()
pil_msk.plot()
#plt.imshow(masked_thinned_MPIL, 'spring', interpolation='none', alpha=1) # 'spring' represents the bright pink color
plt.xlabel('Carrington Longitude [deg]',fontsize = 12)
plt.ylabel('Latitude [deg]',fontsize = 12)
plt.show()
#----

# fig = plt.figure(figsize=(8,6))
# hmi_magmap.plot()
# plt.imshow(masked_thinned_MPIL, 'spring', interpolation='none', alpha=1) # 'spring' represents the bright pink color
# plt.xlabel('Carrington Longitude [deg]',fontsize = 12)
# plt.ylabel('Latitude [deg]',fontsize = 12)
# plt.show()