import drms
import pandas as pd

import os
import numpy as np


client = drms.Client()  
client = drms.Client(email='adithya1@atulbhats.com') # #adithya1@atulbhats.com
query_time_range='hmi.sharp_cea_720s[12160][2024.11.01_12:30:00_TAI-2024.11.01_15:00:00_TAI]' #time range
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

