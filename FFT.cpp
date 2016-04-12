#include<iostream>
#define _USE_MATH_DEFINES
#include<math.h>
#include "FFT.h"
using namespace std;

#define T 50
#define N 60
#define G 128

void FFT::changeOrder(double (**rhok),double (**rhoki),int n,int t)
{
	//Define the variable.
    double temp;
    int k;
    int lh=n/2;
    int j=lh;
    int n1=n-2;

	//Change the order.
    for(int i=1;i<=n1;++i)
    {
         if(i<j)
         {
              temp=rhok[t][i];
              rhok[t][i]=rhok[t][j];
              rhok[t][j]=temp;
              temp=rhoki[t][i];
              rhoki[t][i]=rhoki[t][j];
              rhoki[t][j]=temp;
         }
         k=lh;
         while(j>=k)
         {
              j=j-k;
              k=(int)(k/2+0.5);
         }
         j=j+k;
    }
}

void FFT::FFT_1D(double (**rhoj),double (**rhoji),double (**rhok),double (**rhoki),int len,int t)
{ 
	//Define the variable.
    int m=ceil(log((double)len)/log(2.0));
    int l,b,j,p,k,z;
    double rkb,ikb;
    int n=1<<m;

	//Define the pointer.
    double *rcos;
    double *isin;

	//The memory allocation of the pointer to pointer on the stack.
	rcos=new double[n/2];
	isin=new double[n/2];

	//Assignment.
	z=t;
    for(l=0;l!=n/2;++l)           
    {
         rcos[l]=cos(l*(M_PI)*2/n);
         isin[l]=sin(l*(M_PI)*2/n);
    }

    //memcpy(rhok,rhoj,sizeof(double)*len);
    //memcpy(rhoki,rhoji,sizeof(double)*len);
	for(int i=0;i!=T;++i)
	for(int j=0;j!=G;++j)
	{
		rhok[i][j]=rhoj[i][j];
		rhoki[i][j]=rhoji[i][j];
	}

    for(l=len;l!=n;++l)
    {
         rhok[t][l]=0;
         rhoki[t][l]=0;
    }

	//Change the order
    changeOrder(rhok,rhoki,n,z);

    for(l=1;l<=m;++l)
    {
         b=(int)(pow((double)2,l-1)+0.5);
         for(j=0;j!=b;++j)
         {
              p=j*(int)(pow((double)2,m-l)+0.5);
              for(k=j;k<n;k+=(int)(pow((double)2,l)+0.5))
              {
                   rkb=rhok[t][k+b]*rcos[p]+rhoki[t][k+b]*isin[p];
                   ikb=rhoki[t][k+b]*rcos[p]-rhok[t][k+b]*isin[p];
                   rhok[t][k+b]=rhok[t][k]-rkb;
                   rhoki[t][k+b]=rhoki[t][k]-ikb;
                   rhok[t][k]=rhok[t][k]+rkb;
                   rhoki[t][k]=rhok[t][k]+ikb;
              }
          }
    }

	//After the use of the pointer, release the memory on the stack.
    delete []rcos;
    delete []isin;
}

void FFT::IFFT_1D(double (**phik),double (**phiki),double (**phij),double (**phiji),int len,int t)
{
	//Define the variable
    int m=ceil(log((double)len)/log(2.0));
    double n=1<<m;

    for(int i=0;i!=n;++i)
    phiki[t][i]=-phiki[t][i];

	//FFT_1D() subroutine.
    FFT_1D(phik,phiki,phij,phiji,n,t);

    for(int i=0;i!=len;++i)
    {
         phij[t][i]=phij[t][i]/n;
         phiji[t][i]=-phiji[t][i]/n;
    }
}

