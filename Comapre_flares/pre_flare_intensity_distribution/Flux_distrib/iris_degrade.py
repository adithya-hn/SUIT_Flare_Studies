import matplotlib.pyplot as plt
import numpy as np

y=[]
for i in range(7):
    if i ==0:
        m=1.87e3
        c=3.39e6
    if i ==1:
        m=1.87e3
        c=3.39e6
    if i ==2:
        m=2.32e3
        c=-2.77e6
    if i ==3:
        m=2.84e3
        c=-4.78e6
    if i ==4:
        m=2.84e3
        c=-4.78e6
    if i ==5:
        m=3.28e3
        c=-4.43e6
    if i ==6:
        m=3.28e3
        c=-4.43e6
    k=10000
    y.append(k*m+c)
d=np.array([0.95829199, 0.95652939, 1., 0.92546307, 0.92380995,0.88685412, 0.88459147])
y=np.array(y)
x=np.arange(1,8)
plt.plot(x,y)
plt.plot(x,2e7/d)
plt.show()