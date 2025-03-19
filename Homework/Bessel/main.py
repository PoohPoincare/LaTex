import numpy as np
import matplotlib.pyplot as plt

def BesselJ(x,n,M):
    e=0.05
    J=np.zeros((M),float)
    J[M-1]=0
    J[M-2]=1
    sum=0
    for i in range(M-2,0,-1):
        J[i-1]=(2*i/x)*J[i]-J[i+1]
        if J[i-1]**2>=1/e**2:
            for j in range(i-1,M):
                J[j]=J[j]*e
    sum+=J[0]**2
    for i in range(1,M):
        sum+=2*J[i]**2
    for i in range(0,M):
        J[i]=J[i]/np.sqrt(sum)
    return J[n]
def FunX():
    xx=np.arange(0.1,10,0.1)
    aa,bb,cc,dd,ee,ff=[],[],[],[],[],[]
    for i in xx:
        aa.append(BesselJ(i,0,100))
        bb.append(BesselJ(i,1,100))
        cc.append(BesselJ(i,2,100))
        dd.append(BesselJ(i,3,100))
        ee.append(BesselJ(i,4,100))
        ff.append(BesselJ(i,5,100))
    fig,ax0=plt.subplots()
    ax0.plot(xx,aa,xx,bb,xx,cc,xx,dd,xx,ee,xx,ff)
    ax0.legend(('n=0','n=1','n=2','n=3','n=4','n=5'))
    plt.show()
def FunN():
    nn=np.arange(0,50)
    yy,zz=[],[]
    for k in nn:
        yy.append(BesselJ(1,k,100))
        zz.append(BesselJ(10,k,100))
    fig,ax=plt.subplots()
    ax.plot(nn,yy,nn,zz)
    ax.legend(('x=1','x=10'))
    plt.show()
if __name__=='__main__':
    FunX()
    FunN()

