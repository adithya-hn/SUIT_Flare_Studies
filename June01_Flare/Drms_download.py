import drms
import os
import numpy as np


client = drms.Client()  
client = drms.Client(email='physics@atulbhats.com') # #adithya1@atulbhats.com
query_time_range='hmi.M_45s[2024.11.01_11:00:00_TAI-2024.11.01_15:00:00_TAI]' #time range
keywords = "T_REC,CRPIX1,CRPIX2,CRVAL1,CRVAL2,CDELT1,CDELT2,CUNIT1,CUNIT2,CTYPE1,CTYPE2,DATE-OBS,TIMESYS,TELESCOP,INSTRUME,BUNIT,RSUN_OBS,DSUN_OBS,CRLN_OBS,CRLT_OBS,CAR_ROT,CROTA2" #required headers for sunpy map
query = client.query(query_time_range,key=keywords) # just for query
print(query)

export_request = client.export(query_time_range) # No proper headers will be there but fast download
#export_request = client.export(query_time_range, method='url', protocol='fits') #works perfectly if jsoc server back in action, header will be proper

# Save the result to a CSV file
csv_file = f'drms_query_results_{query_time_range}.csv' #header information for updating it later
query.to_csv(csv_file, index=False)
print(f"Query results saved to {csv_file}")

out_dir='/Analysis/Projects_Data/Flare_Data/Nov01_Flare_Data2/HMI'
export_request.download(out_dir,verbose=True)

