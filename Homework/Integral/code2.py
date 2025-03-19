import numpy as np
def H(x,y):
    if x**2+y**2>=1:
        return 0
    else:
        return 1
def MonteCarlo(N):
    j=0
    x=-1+2*np.random.rand(N)
    y= -1+2*np.random.rand(N)
    z=np.random.rand(N)
    for i in range(0,N):
        if H(x[i],y[i])>=z[i]:
            j+=1
    I=j/N*4
    return I