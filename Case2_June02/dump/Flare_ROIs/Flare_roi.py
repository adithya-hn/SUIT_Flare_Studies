import suitkit
import numpy as np
from datetime import datetime

csv_file = '/data/sreejith/MCNS_POC/Flare_ROIs/Database/SUIT_Data_Jun_catalogue.csv'
start_time = datetime(2024, 6, 2, 2, 50, 0)
end_time = datetime(2024, 6, 2, 5, 10, 0)
filtered_files = suitkit.list_files_in_time_range(csv_file, start_time, end_time,img_type='ROI')
np.savetxt('Flare_files_Jun2_M1.2.dat',filtered_files , delimiter=',', fmt='%s')

suitkit.map_movie(filtered_files,Filter='NB03',start_time=start_time,end_time=end_time)