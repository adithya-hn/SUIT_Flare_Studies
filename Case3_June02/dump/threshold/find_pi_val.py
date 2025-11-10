import numpy as np

import matplotlib.pyplot as plt
a=2
circ_pts=[]
sqr_ptx=[]
for i in range(10):
    x=np.random.uniform(-a,a,1)[0]
    y=np.random.uniform(-a,a,1)[0]
    print(x)
    if x**2 +y**2 <= a:
        circ_pts.append((x,y))
    else:
        sqr_ptx.append((x,y))





