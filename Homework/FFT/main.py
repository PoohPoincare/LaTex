import numpy  as np
import matplotlib.pyplot as plt

def DataRev(data):
    N=len(data)
    data_rev=[]
    width=int(np.log(N)/np.log(2))#计算位宽 width，即 nmax 表示的二进制数的位数。
    for i in range(0,N):
        j= '{:0{width}b}'.format(i, width=width)#使用格式化字符串 '{:0{width}b}' 将整数 i 转换为 width 位的二进制字符串，前面用 0 补齐。
        data_rev.append(data[int(j[::-1], 2)])#将生成的二进制字符串 j 反转（j[::-1]），然后将其转换回整数并返回。
    return data_rev
def Coeff(N,k,isign):
    return  np.exp(isign*2*np.pi*1j*k/N)
def FFT(data,isign):
    data=np.array(data)
    data=DataRev(data)
    Length=len(data)
    Width=int(np.log(Length)/np.log(2))
    for i in range(Width):
        for j in range(0,Length,2**(i+1)):
            for k in range(0,2**i):
                x=k+j
                y=x+2**i
                W = Coeff(2**(i+1),k,isign)
                data[x],data[y]= data[x]+W*data[y], data[x]-W*data[y]
    if isign==-1:
        i=0
        while(i<Length):
            data[i]=data[i]/Length
            i+=1
    return data
def bessl(z,M):
    N=2**M
    data=[]
    for i in range(0,N):
        data.append(np.cos(z*np.cos(2*np.pi*i/N))+np.sin(z*np.cos(2*np.pi*i/N))*1j)
        #data.append(np.exp(z*np.cos(2*np.pi*i/N)))
    data=FFT(data,1)
    for i in range(0,N):
        data[i]=data[i]/(N*(1j)**i)
        data[i]=data[i].real
    return data
if __name__=='__main__':
    J1=bessl(1,12)
    j1=[]
    xx=np.arange(0,50)
    for i in range(0,50):
      j1.append(J1[i])
    plt.plot(xx,j1)
    plt.show()

