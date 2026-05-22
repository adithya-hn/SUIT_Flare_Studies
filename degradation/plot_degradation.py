import numpy as np
import matplotlib.pyplot as plt



data=(np.loadtxt('nb3_degradation.csv',delimiter=',')).transpose( )
dt=np.array(data[0],dtype=float)
vals=np.array(data[1],dtype=float)
plt.scatter(dt,vals)
plt.show()