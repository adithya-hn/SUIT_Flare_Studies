import numpy as np
import pandas as pd


Helios=(np.load("cdte_data_flare_7.npy", allow_pickle=True)).transpose()
er=np.sqrt(np.array(Helios[3]**2, dtype=np.float64)+np.array(Helios[4]**2, dtype=np.float64))
np.savetxt('helios_CdTe_c7.csv',np.c_[Helios[0],(Helios[1]+Helios[2]),er],delimiter=',',header='Time,Total,CdTe1+2Er',comments='',fmt='%s')
