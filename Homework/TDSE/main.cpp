#include<iostream>
#include<cmath>
#include<vector>
#include <bitset>
#include <string>//to_string()
#include <algorithm>//reverse()
#include <fstream>
#define PI 3.14159
#define NN 32768
#define Wi 15
using namespace std;
class Wave
{
    public:
        float sigma,delta,k;
        float x0=0,y0=0;
        int istart=0;
        int Nx,Ny,Nh;
        int length=5,width=5;
        //length和width分别是势垒的厚度和宽度的一半，单位是delta。需要注意，取值应该是python版本的一半。
        float Vmax,T;
        float *X;
        float *Y;
        float *potential;
        float **psi;


        int initialize(float x_min,float x_max,float y_min,float y_max,float in_delta,float k,float Vmax,float Ek,int in_istart);
        int operate(float vector[]);
        int TimeInt(float An_r[],float An_i[],int m);
    
};
int Wave::initialize(float x_min,float x_max,float y_min,float y_max,float in_delta,float in_k,float in_Vmax,float Ek,int in_istart)
{   T=Ek;
    istart=in_istart;
    Vmax=in_Vmax;
    delta=in_delta;
    sigma=20*delta;
    k=in_k;
    Nx=int((x_max-x_min)/delta);
    Ny=int((y_max-y_min)/delta);
    cout<<"grid number"<<Nx<<'|'<<Ny<<endl;
    Nh=Nx*Ny;
    X =new float[Nx];
    Y =new float[Ny];
    for(int i=0;i<Nx;i++)
    {
        X[i]=x_min+i*delta;
        //cout<<X[i]<<endl;
    }
    for(int j=0;j<Ny;j++)
    {
        Y[j]=y_min+j*delta;
    }
    psi=new float*[2];
    psi[0]=new float[Nh]();
    psi[1]=new float[Nh]();
    if(istart==0){
        for(int i=0;i<Ny;i++)
    {
        for(int j=0;j<Nx;j++)
        {
            psi[0][i*Nx+j]=exp(-pow((X[j]-x0),2)/(4*pow(sigma,2))-pow((Y[i]-y0),2)/(4*pow(sigma,2)))*cos(k*X[j]);
            //psi[0][i*Nx+j]=1;
           //cout<<exp(-pow(0,2)/(4*pow(sigma,2))-pow(0,2)/(4*pow(sigma,2)))<<endl;
            psi[1][i*Nx+j]=exp(-pow((X[j]-x0),2)/(4*pow(sigma,2))-pow((Y[i]-y0),2)/(4*pow(sigma,2)))*sin(k*X[j]);            
        }
    }
    }
    else{
        ifstream infile; 
        infile.open("Psi_r.dat");
        for(int i=0;i<Nh;i++)
        {   
        infile>>psi[0][i];
        }
        infile.close();
        infile.open("Psi_i.dat");
        for(int i=0;i<Nh;i++)
        {   
        infile>>psi[1][i];
        }
        infile.close();        
    }

    //cout<<X[100]<<endl;
    //cout<<exp(-pow((X[100]-x0),2)/(4*pow(sigma,2))-pow((Y[100]-y0),2)/(4*pow(sigma,2)))*cos(k*X[100])<<endl;
    potential= new float[Nh]();
    int x1=Nx/2-length;
    int x2=Nx/2+length;
    int y1=Ny/2-width;
    int y2=Ny/2+width;
    for(int j=x1;j<x2;j++)
    {
        for(int i=0;i<Ny;i++)
        {
            potential[i*Nx+j]=Vmax;
        }        
        for(int i=y1;i<y2;i++)
        {
            potential[i*Nx+j]=0;
        }
    }
    delete[] X;
    delete[] Y;
    return 0;
}
int Wave::operate(float vector[])
{   
    float *Temp=new float[Nh]();
    Temp[0]=(potential[0]+4*T)*vector[0]-T*(vector[1]+vector[Nx]);
    for(int i=1;i<Nx;i++)
    {
        Temp[i]=(potential[i]+4*T)*vector[i]-T*(vector[i-1]+vector[i+1]+vector[i+Nx]);
    }
    for(int i=Nx;i<Nh-Nx;i++)
    {
        Temp[i]=(potential[i]+4*T)*vector[i]-T*(vector[i-1]+vector[i+1]+vector[i+Nx]+vector[i-Nx]);
    }
    for(int i=Nh-Nx;i<Nh-1;i++)
    {
        Temp[i]=(potential[i]+4*T)*vector[i]-T*(vector[i-1]+vector[i+1]+vector[i-Nx]);
    }
    Temp[Nh-1]=(potential[Nh-1]+4*T)*vector[Nh-1]-T*(vector[Nh-2]+vector[Nh-1-Nx]);
    for(int i=0;i<Nh;i++)
    {
        vector[i]=Temp[i];
    }
    delete[] Temp;
    return 0;
}
int Wave::TimeInt(float An_r[],float An_i[],int m)
{
    float **Temp=new float*[2];
    Temp[0]=new float[Nh](),Temp[1]=new float[Nh]();
    float **Temp0=new float*[2];
    Temp0[0]=new float[Nh](),Temp0[1]=new float[Nh]();
    float **Temp1=new float*[2];
    Temp1[0]=new float[Nh](),Temp1[1]=new float[Nh]();
    for(int i=0;i<Nh;i++)
    {
        Temp0[0][i]=psi[0][i];
        Temp0[1][i]=psi[1][i];
    }
    operate(psi[0]);
    operate(psi[1]);
    for(int i=0;i<Nh;i++)
    {
        Temp1[0][i]=psi[0][i];
        Temp1[1][i]=psi[1][i];
        psi[0][i]=0.5*(An_r[0]*Temp0[0][i]-An_i[0]*Temp0[1][i])+(An_r[1]*Temp1[0][i]-An_i[1]*Temp1[1][i]);
        psi[1][i]=0.5*(An_i[0]*Temp0[0][i]+An_r[0]*Temp0[1][i])+(An_r[1]*Temp1[1][i]+An_i[1]*Temp1[0][i]);
    }
    float *temp_r=new float[Nh];
    float *temp_i=new float[Nh];
    for(int j=2;j<m;j++)
    {
        for(int i=0;i<Nh;i++)
        {
            temp_r[i]=Temp1[0][i];
            temp_i[i]=Temp1[1][i];
        }
        operate(Temp1[0]);
        operate(Temp1[1]);   
        for(int i=0;i<Nh;i++)
        {
            Temp[0][i]=2*Temp1[0][i]-Temp0[0][i];
            Temp[1][i]=2*Temp1[1][i]-Temp0[1][i];
            psi[0][i]=psi[0][i]+(An_r[j]*Temp[0][i]-An_i[j]*Temp[1][i]);
            psi[1][i]=psi[1][i]+(An_r[j]*Temp[1][i]+An_i[j]*Temp[0][i]);
            Temp1[0][i]=Temp[0][i];
            Temp1[1][i]=Temp[1][i];
            Temp0[0][i]=temp_r[i];
            Temp0[1][i]=temp_i[i];
        }
        if(j%200==0)
        {   
            float fm=float(m);
            cout<<'#';
        }

    }
    delete [] Temp0;
    delete[] Temp1;
    delete[] Temp;
    return 0;

}
float DataRev(float data[])
{   
    float data_rev[NN];
    int width=int(log(NN)/log(2));
    for(unsigned int i=0;i<NN;i++)
    {   
        bitset<Wi> bit(i);//i转二进制数
        string iStr=bit.to_string();//二进制数转字符串
        reverse(iStr.begin(), iStr.end());//字符串倒序
        int j = stoi(iStr, nullptr, 2);//字符串转十进制数
        data_rev[i]=data[j];
        //cout<<iStr<<endl;
    }
    for(int i=0;i<NN;i++)
    {
        data[i]=data_rev[i];
    }
    return 0;
}
float Coeff_r(int N,int k,int isign)
{
    return cos(isign*2*3.1415926*k/N);
}
float Coeff_i(int N,int k,int isign)
{
    return sin(isign*2*3.1415926*k/N);
}
int FFT(float data_r[],float data_i[], int isign,int N,int Width)
{   cout<<"begin DataRev"<<endl;
    DataRev(data_r);
    DataRev(data_i);
    for(int i=0;i<Width;i++)
    {
        for(int j=0;j<N;j+=pow(2,i+1))
        {
            for(int k=0;k<pow(2,i);k++)
            {
                int x=k+j;
                int y=x+pow(2,i);
                float W_r = Coeff_r(pow(2,i+1),k,isign),W_i = Coeff_i(pow(2,i+1),k,isign);
                float temp_r=data_r[y],temp_i=data_i[y];
                data_r[y]=data_r[x]-W_r*temp_r+W_i*temp_i;
                data_i[y]=data_i[x]-W_r*temp_i-W_i*temp_r;
                data_r[x]=data_r[x]+W_r*temp_r-W_i*temp_i;
                data_i[x]=data_i[x]+W_r*temp_i+W_i*temp_r;
            }

        }
    }
    if(isign==1)
    {
        for(int i=0;i<N;i++)
        {
            data_r[i]=2*data_r[i]/N;
            data_i[i]=2*data_i[i]/N;
        }
    }
    return 0;
}
int Bessel(float data_r[],float data_i[],float z,int N,int M)
{   
    for(int i=0;i<N;i++)
    {
        data_r[i]=cos(-z*cos(2*PI*i/N));
        data_i[i]=sin(-z*cos(2*PI*i/N));
    }
    cout<<"begin FFT"<<endl;
    FFT(data_r,data_i,1,N,M);
    cout<<"FFT down"<<endl;
    /*
    for(int i=0;i<N;i=i+4)

    {   data_r[i]=2*data_i[i];
        data_r[i+1]=2*data_i[i+1];
        data_r[i+2]=-2*data_r[i+2];
        data_r[i+3]=-2*data_i[i+3];
    }
    */
    return 0;
}

int main()
{   
    cout<<"initialize:...";
    Wave Psi;
    Psi.initialize(-10,20,-10,10,0.1,5,0.1,0.1,0);
    cout<<"down"<<endl;
    cout<<"computing coeffient an:...";
    float An_r[NN];
    float An_i[NN];
    int m=5000,timesteps=1100;
    Bessel(An_r,An_i,timesteps,NN,Wi);
    cout<<"down"<<endl;
    cout<<"begin TimeInt:";
    //time integral
    Psi.TimeInt(An_r,An_i,m);
    cout<<endl<<"down"<<endl;
    cout<<"write dat"<<endl;
    ofstream outfile;
    outfile.open("Psi_r.dat");
    for(int i=0;i<Psi.Nh;i++)
    {   

        outfile<<Psi.psi[0][i]<<endl;
    }
    outfile.close();
    outfile.open("Psi_i.dat");
    for(int i=0;i<Psi.Nh;i++)
    {   
        outfile<<Psi.psi[1][i]<<endl;
    }
    outfile.close();
    cout<<"completed"<<endl;
    
}