from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from sunpy.time import parse_time
from sys import path as sys_path
sys_path.append('/home/adithya/Adithya_repos')
from plots_styl import set_pub_style
set_pub_style()
from scipy.optimize import curve_fit


data1=np.loadtxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/degradation/qual_data_NB04_bk1.csv',delimiter=',',dtype='str',skiprows=1)
data2=np.loadtxt('/Analysis/Research_Projects/Flare_studies/SUIT_Flares/degradation/qual_data_NB04_bk2.csv',delimiter=',',dtype='str',skiprows=1)

dt1_ =np.array(data1[:,0],dtype='datetime64')
idx1 = np.argsort(dt1_)
dt1 = dt1_[idx1]
val1=np.array(data1[:,1][idx1],dtype=float)
t_sec1=(np.array(data1[:,0],dtype='datetime64[s]').astype('timedelta64[s]').astype(float))


dt2_ =np.array(data2[:,0],dtype='datetime64')
idx2 = np.argsort(dt2_)
dt2 = dt2_[idx2]
val2=np.array(data2[:,1][idx2],dtype=float)
t_sec2=(np.array(data2[:,0],dtype='datetime64[s]').astype('timedelta64[s]').astype(float))



def linear(x, a, b):
    return a*x + b

def expo(x, A, B):
    return A * np.exp(B * x)


p_lin1, _ = curve_fit(linear, t_sec1, val1)
y_lin1 = linear(t_sec1, *p_lin1)

p_lin2, _ = curve_fit(linear, t_sec2, val2)
y_lin2 = linear(t_sec2, *p_lin2)


'''
plt.figure(figsize=(14,8))
plt.title('NB04 degradation')
plt.scatter(dt,np.array(data[:,1][idx],dtype=float),marker='+',color='r')
plt.xlabel('Time (UT)')
plt.ylabel('Exposure normalized relative intensity')
#plt.plot(dt,y_exp,label='exponential')
plt.plot(dt,y_lin,label='linear')
plt.plot(dt, np.sum((val - y_fit)**2))
#time_formatter = mdates.DateFormatter('%Y')  # Format as HH:MM
#plt.gca().xaxis.set_major_formatter(time_formatter)
plt.savefig('baking_1_degradation.png',dpi=300)
plt.show()

'''

residuals1=(val1 - y_lin1)
residuals2=(val2 - y_lin2)
# --- figure -------------------------------------------------


#gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.05)
fig, ax1 = plt.subplots(2, 1,height_ratios=[3, 1], sharex=True, figsize=(18,12),gridspec_kw={'hspace': 0})  # no vertical spacing 

# --- top panel: fit ----------------------------------------
ax1[0].scatter(dt1,np.array(data1[:,1][idx1],dtype=float),marker='+',color='r')
ax1[0].plot(dt1,y_lin1,label=f'linear fit (slope ({p_lin1[0]:.3}))')

ax1[0].scatter(dt2,np.array(data2[:,1][idx2],dtype=float),marker='+',color='r')
ax1[0].plot(dt2,y_lin2,label=f'linear fit (slope ({p_lin2[0]:.3}))')

ax1[0].legend()

ax1[1].set_xlabel('Time (UT)')
ax1[0].set_ylabel('Exposure normalized relative intensity')

# --- bottom panel: residuals -------------------------------
ax1[1].axhline(0, color='k', linewidth=0.8)
ax1[1].plot(dt1, residuals1, 'o')
ax1[1].plot(dt2, residuals2, 'o')
ax1[1].set_ylabel('Residuals')
#ax1.xaxis.get_offset_text().set_visible(False)


# remove x tick labels from top panel
plt.setp(ax1[0].get_xticklabels(), visible=False)

#plt.tight_layout()
plt.savefig('plot.png',dpi=300)
plt.close()

obs_time_stps1=np.array(["2024-07-10T03:59:00","2024-07-10T13:37:00"],dtype='datetime64') #"2024-06-02T02:50:00","2024-06-02T06:50:00",
obs_time_stps2=np.array(["2024-10-08T23:56:00","2024-11-01T00:16:00","2024-11-01T12:31:00","2024-11-12T22:22:00","2024-11-13T15:08:00"],dtype='datetime64')


obs_dt_sec1=obs_time_stps1.astype('timedelta64[s]').astype(float)
obs_dt_sec2=obs_time_stps2.astype('timedelta64[s]').astype(float)

y_lin2 = np.array(linear(obs_dt_sec2, *p_lin2))
y_lin1 = np.array(linear(obs_dt_sec1, *p_lin1))


th_put=np.concatenate((y_lin1, y_lin2))
print(th_put)

norm_th=th_put/th_put.max()
print(norm_th)

qs_median=np.array([5203,5276,5661,4764,4793,4373,4355])
norm_qsmedian=qs_median/qs_median.max()
print(norm_qsmedian)