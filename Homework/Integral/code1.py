import numpy as np
def f(x):
    return 4/(1+x**2)
def trapzoidal(N):
    x1,x2=0.,1.
    I=0
    d=(x2-x1)/N
    x=np.arange(x1,x2+d,d)
    for i in range(0,N):
        m=(x[i]+x[i+1])/2
        I+=f(m)*d
    return I
def simpson(N):
    x1,x2=0.,1. 
    I=0
    d=(x2-x1)/N
    x=np.arange(x1,x2+d,d)
    for i in range(0,N):
        m=(x[i]+x[i+1])/2
        I+=(d/6)*(f(x[i])+f(m)*4+f(x[i+1]))
    return I