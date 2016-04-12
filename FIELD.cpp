#include<iostream>
#include "FIELD.h"
#include "FFT.h"
#define _USE_MATH_DEFINES
#include<math.h>
using namespace std;

#define T 50
#define N 5000
#define G 128

void FIELD::field(double *Xj,double (**rhoj),double (**rhoji),double (**ej),double (**eji),double epsi,double dx,int t)
{
    //Define the variable.
    FFT fft;
    int ng=128;
    int ng1=64;
    int ng2=128;
    double kdx2[128];
    double sm[128];
    double ksqi[128];
    //Define the pointer to pointer.
    double **rhok;
    double **rhoki;
    double **phik;
    double **phiki;
    double **phij;
    double **phiji;
    int v;
    
    //Memory allocation of the pointer to pointer on the stack.
    rhok=new double*[T];
    rhoki=new double*[T];
    phik=new double*[T];
    phiki=new double*[T];
    phij=new double*[T];
    phiji=new double*[T];
    
    for(int i=0;i!=T;++i)
    {
    	rhok[i]=new double[G];
    	rhoki[i]=new double[G];
    	phik[i]=new double[G];
    	phiki[i]=new double[G];
    	phij[i]=new double[G];
    	phiji[i]=new double[G];
    }
    
    //Assignment
    v=t;
    for(int k=0;k!=G;++k)
    {
    	kdx2[k]=(2*M_PI/ng)*k+0.0005;
    	sm[k]=0;
    	ksqi[k]=0;
    }
    //Caculate the smooting factor sm and 1/(K*K).
    for(int k=0;k!=G;++k)
    {
    	sm[k]=exp(2*sin(kdx2[k])*sin(kdx2[k])-2*tan(kdx2[k])*tan(kdx2[k])*tan(kdx2[k])*tan(kdx2[k]));
    	ksqi[k]=epsi*sm[k]*sm[k]/((2*sin(kdx2[k])/dx)*(2*sin(kdx2[k])/dx));
    }
    
    for(int t=0;t!=T;++t)
    	for(int k=0;k!=G;++k)
    	{
    		rhok[t][k]=0.0;
    		rhoki[t][k]=0.0;
    		phik[t][k]=0.0;
    		phiki[t][k]=0.0;
    		phij[t][k]=0.0;
    		phiji[t][k]=0.0;
    	}
    
    //Caculate the rhok[50][128] via rhoj[50][128].
    //FFT_1D() subroutine
    fft.FFT_1D(rhoj,rhoji,rhok,rhoki,G,v);
    
    //Caculate the phik[t][128] and phiki[t][128] via rhok[t][128] and rhoki[t][128].
    for(int k=0;k!=G;++k)
    {
    phik[v][k]=rhok[v][k]*ksqi[k];
    phiki[v][k]=rhoki[v][k]*ksqi[k];
    }
    
    //IFFT_1D subroutine
    fft.IFFT_1D(phik,phiki,phij,phiji,G,v);
    
    //Caculate the ej[t][128] via phij[t][128].
    for(int j=0;j!=G;++j)
    {
    	if(j==0)
    		ej[v][j]=(phij[v][j+127]-phij[v][j+1])/(2*dx);
    	else if(j==127)
    		ej[v][j]=(phij[v][j-1]-phij[v][j-127])/(2*dx);
    	else
    	    ej[v][j]=(phij[v][j-1]+phij[v][j+1])/(2*dx);
    }
    
    //After the use of the pointer, release the memory on the stack.
    for(int t=0;t!=T;++t)
    {
    	delete []rhok[t];
    	delete []rhoki[t];
    	delete []phik[t];
    	delete []phiki[t];
    	delete []phij[t];
    	delete []phiji[t];
    }
    delete []rhok;
    delete []rhoki;
    delete []phik;
    delete []phiki;
    delete []phij;
    delete []phiji;
}
