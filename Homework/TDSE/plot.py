import numpy  as np
from sys import version
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import os
def grid(min,max,delta):
    #生成网格,从min到max，间距为delta
    line=np.arange(min,max,delta)
    return line
def gragh(X,Y,V):
    #绘制波函数的3D图像
    fig = plt.figure()
    ax = Axes3D(fig)
    fig.add_axes(ax)
    print(len(X),len(Y),len(V))
    Z=np.zeros((len(Y),len(X)),float)
    for i in range(0,len(Y)):
        for j in range(0,len(X)):
            Z[i][j]=V[i*len(X)+j]
    X, Y = np.meshgrid(X, Y)
    ax.plot_surface(X, Y, Z,rstride=1, cstride=1, cmap=cm.viridis)
    plt.show()
def level(X,Y,V,psi):
#波函数的俯视图
    plt.style.use('_mpl-gallery-nogrid')
    Z1=np.zeros((len(Y),len(X)),float)
    Z2=np.zeros((len(Y),len(X)),float)
    #y在前，x在后的二维数组
    for i in range(0,len(Y)):
        for j in range(0,len(X)):
            Z1[i][j]=psi[i*len(X)+j]
            Z2[i][j]=V[i*len(X)+j]
    level1 = np.linspace(Z1.min(), Z1.max(),10)
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.contourf(X, Y, Z1,  levels=level1,cmap=cm.turbo)
    ax.contourf(X, Y, Z2,  levels=level1,cmap=cm.turbo)
    plt.show()
def potential(X,Y,Vmax,length,width,d):
    #生成一个中间带有小孔的势垒
    Ny=len(Y)
    Nx=len(X)
    Nh=Nx*Ny
    V=np.zeros((Nh),float)
    mx=int(len(X)/2)
    my=int(len(Y)/2)
    #mx,my是整个网格的中间位置，length和width对应小孔在x和y方向的尺度
    x1=int((X[mx]-X[0]-length/2)/d)
    x2=int((X[mx]-X[0]+length/2)/d)
    y1=int((Y[my]-Y[0]-width/2)/d)
    y2=int((Y[my]-Y[0]+width/2)/d)
    for j in range(x1,x2):
        for i in range(0,Ny):
            V[i*Nx+j]=Vmax
        for i in range(y1,y2):
            V[i*Nx+j]=0.
    return V
if __name__=='__main__':
    '''
    data_r = pd.read_csv('TDSE/psi_r.dat', sep='\t', dtype=float)
    data_i = pd.read_csv('TDSE/psi_i.dat', sep='\t', dtype=float)
    data_r = np.fromfile('TDSE/psi_r.dat', dtype=np.float32)
    data_i = np.fromfile('TDSE/psi_i.dat', dtype=np.float32)

    '''
    path=os.path.join(os.getcwd(),'TDSE/psi_r.dat')
    f=open(path)
    print("read file")
    line=f.readline().strip()
    data_r=[]
    print(len(data_r))
    i=0
    while line:
        data_r.append(line)
        line=f.readline().strip()
        i+=1
    print(i)
    f.close()
    f=open('TDSE/psi_i.dat')
    print("read file")
    line=f.readline().strip()
    data_i=[]
    i=0
    while line:
        data_i.append(line)
        line=f.readline().strip()
        i+=1
    print(i)
    f.close()
    
    psi_r=[]
    psi_i=[]
    psi2=[]
    print(len(data_i))
    for j in range(0,len(data_i)):
           psi_r.append(float(data_r[j]))
           psi_i.append(float(data_i[j]))
           psi2.append(psi_r[j]**2+psi_i[j]**2)
    d=0.1
    X=grid(-10,20,d)
    Y=grid(-10,10,d)
    print("prepare figure")
    V=potential(X,Y,0.5,10*d,10*d,d)
    level(X,Y,V,psi2)
    #gragh(X,Y,psi2)

