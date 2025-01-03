from astropy.io import fits
import sunpy.map
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from datetime import datetime

hlios=fits.open('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/June01_Flare/Helios/02/STH_2024_06_02.fits')
header=hlios[1].header
time_counter=header['TTYPE3']
#print(time_counter)
Data=hlios[1].data

#print(Data)
Data_=Data.transpose()
x1=Data[10][1]
x2=Data[11][1]
x3=Data[12][1]
print(x1,x2,x3)
cdt1=[]
cdt2=[]
date_array=[]

for i in range (len(Data)):
    cdt1.append(Data_[i][3])
    cdt2.append(Data_[i][4])

    #date_array.append(Time(Data_[i][1], scale='utc', format='isot') )
    date_array.append(datetime.strptime((Data_[i][1])[:26], '%Y-%m-%dT%H:%M:%S.%f'))





print(len(date_array),len(cdt1))
plt.plot(date_array,cdt1)
plt.plot(date_array,cdt2)
plt.yscale('log')
plt.show()

binned_counts1=[]
binned_counts2=[]
binned_date=[]
for i in range(0, len(cdt1), 60):
    # Take the next 10 elements from the list of counts
    count_batch1 = cdt1[i:i+60]
    count_batch2 = cdt2[i:i+60]

    # Average the counts in the batch
    binned_counts1.append(np.sum(count_batch1))
    binned_counts2.append(np.sum(count_batch2))

    # Take the 10th timestamp (if available)
    if i + 59 < len(date_array):
        binned_date.append(date_array[i + 59])

binned_date.insert(0,date_array[0])


print(len(binned_counts1),len(binned_date))
#plt.plot(binned_date,binned_counts1)
plt.plot(binned_date,binned_counts2)
plt.yscale('log')
plt.show()
