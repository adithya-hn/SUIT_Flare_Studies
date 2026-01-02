import drms
import pandas as pd

# c = drms.Client()   # JSOC client

# sharp_series = "hmi.sharp_cea_720s"
# sharp_num = 11487

# series = "hmi.sharp_cea_720s"
# # t = "2024.07.10_04:00:00_TAI"

# # res = c.query(f"{series}[{t}]", key="HARPNUM, NOAA_AR, T_REC")
# # print(res)


# t_start = "2024.07.10_04:00:00_TAI"
# t_end   = "2024.07.10_06:00:00_TAI"

# query_string = f"{sharp_series}[{sharp_num}][{t_start}-{t_end}]"

# # Request only keywords (metadata), not images
# keys = "T_OBS, TOTUSJH, TOTUSJZ, USFLUX, MEANGAM, MEANJZD, AREA_ACR"

# res = c.query(query_string, key=keys)

# # df = pd.DataFrame(res)
# # print(df)

# res = c.query(query_string, key="__all__")
# df = pd.DataFrame(res)
# df.to_csv('SHARP_data.csv',index=False)

import drms
import os
import numpy as np


client = drms.Client()  
client = drms.Client(email='adithya1@atulbhats.com') # #adithya1@atulbhats.com
query_time_range='hmi.sharp_cea_720s[11487][2024.07.10_13:00:00_TAI-2024.07.10_16:14:00_TAI]' #time range
#query_time_range='hmi.M_45s [2024.06.02_06:30:00_TAI-2024.06.02_09:20:00_TAI]' 
keywords = "T_REC,CRPIX1,CRPIX2,USFLUX,MEANSHR,MEANPOT,TOTPOT,SHRGT45,MEANJZD,TOTUSJZ,MEANGBH" #required headers for sunpy map
query = client.query(query_time_range,key=keywords) # just for query
print(query)

export_request = client.export(query_time_range) # No proper headers will be there but fast download
#export_request = client.export(query_time_range, method='url', protocol='fits') #works perfectly if jsoc server back in action, header will be proper

# Save the result to a CSV file
csv_file = f'sharp_headers.csv' #header information for updating it later
query.to_csv(csv_file, index=False)
print(f"Query results saved to {csv_file}")

# out_dir='/Analysis/Research_Projects/Flare_studies/SUIT_Flares/Case4_July10/data/SHARP_DATA'
# export_request.download(out_dir)
