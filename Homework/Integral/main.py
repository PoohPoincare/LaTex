import numpy as np
import matplotlib.pyplot as plt
def f(x):
    return 4/(1+x**2)
def H(x,y):
    if x**2+y**2>=1:
        return 0
    else:
        return 1
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


def MonteCarlo(N):
    j=0
    x=-1+2*np.random.rand(N)
    y= -1+2*np.random.rand(N)
    for i in range(0,N):
        if H(x[i],y[i])>0:
            j+=1
    I=4*j/N
    return I


def err1(N):
    return np.abs(trapzoidal(N)-np.pi)/np.pi
def err2(N):
    return np.abs(simpson(N)-np.pi)/np.pi
def err3(N):
    return np.abs(MonteCarlo(N)-np.pi)/np.pi
def fit():
    fig,ax1=plt.subplots()
    xx=np.arange(0,7,0.5)
    xx1=np.arange(0,9,1)
    y10,y1=[],[]
    for i in range(0,14):
        y10.append(trapzoidal(int(10**(i/2))))
        y1.append(np.log10(err1(int(10**(i/2)))))
    y20,y2=[],[]
    for i in range(0,14):
        y20.append(simpson(int(10**(i/2))))
        y2.append(np.log10(err2(int(10**(i/2)))))
    y30,y3=[],[]
    for i in range(0,9):
        y30.append(MonteCarlo(10**(i)))
        y3.append(np.log10(err3(10**(i))))  
#axisx=np.log10(xx)
    ax1.set_yticklabels(['$10^{-16}$','$10^{-14}$','$10^{-12}$','$10^{-10}$','$10^{-8}$','$10^{-6}$','$10^{-4}$','$10^{-2}$']) 
    ax1.set_xticklabels(['1','0.1','$10^{-2}$','$10^{-3}$','$10^{-4}$','$10^{-5}$','$10^{-6}$','$10^{-7}$'])  
    ax1.plot(xx,y1,'r',xx,y2,'b')
    ax1.plot(xx,y1,'ro',xx,y2,'bo')
    ax1.set_title('relative error of  Trapezoidal rule and  Simpson’s rule')
    plt.legend( ('Trapezoidal',' Simpson’s '),loc='upper right')
    plt.xlabel('grid space')
    plt.ylabel('relative error')
    fig,ax2=plt.subplots() 
    ax2.plot(xx1,y3,'g',xx1,y3,'go')
    ax2.set_yticklabels(['$10^{-4.5}$','$10^{-4}$','$10^{-3.5}$','$10^{-3}$','$10^{-2.5}$','$10^{-2}$', '$10^{-1.5}$','$10^{-1}$','$10^{-0.5}$']) 
    plt.legend( ('Monte-carlo',),loc='upper right')
    plt.xlabel('samples N')
    plt.ylabel('relative error')
    ax2.set_xticklabels(['1','10','$10^{2}$','$10^{3}$','$10^{4}$','$10^{5}$','$10^{6}$','$10^{7}$','$10^{8}    $'])
    #ax2.set_title('relation between the relative error and the MonteCarlo samples N. ')
    fig,ax3=plt.subplots()
    ax3.plot(xx,y10,'r',xx,y20,'b')
    ax3.plot(xx,y10,'ro',xx,y20,'bo')
    fig,ax4=plt.subplots()
    ax4.plot(xx1,y30,'g',xx1,y30,'go')
    plt.show()
if __name__=='__main__':
    fit()
    
    
   

