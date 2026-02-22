

import numpy as np
import matplotlib.pyplot as plt



goes  = (np.loadtxt('goes_xray_lightcurve.csv',delimiter=',',skiprows=1,dtype='str')).transpose() 

goes_dt=np.array(goes[0],dtype='datetime64')
xrs_a=np.array(goes[1],dtype=float)
xrs_b=np.array(goes[2],dtype=float)

plt.plot(goes_dt,xrs_a)
plt.plot(goes_dt,xrs_b)
plt.close()

R=xrs_a/xrs_b

T=  2.7460+ 129.47*R - 966.28*pow(R,2)+ 5517.5*pow(R,3)- 1.8664*pow(10,4)*pow(R,4)+ 3.5951*pow(10,4)*pow(R,5)- 3.6099 *pow(10,4)*pow(R,6)+ 1.4687*pow(10,4)*pow(R,7) 
xrs_b_model=6.9497- 6.0827*T+  1.7364*pow(T,2)- 0.15594*pow(T,3)+  6.7848*pow(10,-3)*pow(T,4)- 1.4446*pow(10,-4)*pow(T,5)* + 1.2089* pow(10,-6)*pow(T,6)
em=xrs_b/xrs_b_model

fig,ax=plt.subplots()
ax2=ax.twinx()
ax2.plot(goes_dt,T,'r',label='Temperature')
#ax.plot(goes_dt,em,label='Emission measure')
ax.set_yscale('log')
plt.grid()
plt.legend()
plt.show()
goes_str=np.array(goes_dt,dtype='str')
np.savetxt('tn_woods_goes_temp_em.csv',np.c_[goes_str,T,em],delimiter=',',header='date_time,temperature,em10^49',fmt='%s',comments='')