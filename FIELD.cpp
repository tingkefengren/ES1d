#include<iostream>
#include<fstream>
#include<stdlib.h>
#include "FIELD.h"
#include "FFT.h"
#define _USE_MATH_DEFINES
#include<math.h>
using namespace std;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

extern const int G;

void FIELD::field(double *Xj,double (**rhoj),double (**rhoji),double (**ej),double epsi,double dx,int t,char *argv1,char *argv2)
{
	double a1=atoi(argv1);			//a1--compensation factor
	double a2=atoi(argv2);			//a2--smoothing factor
	double k[G];				//k--wave number
	double sm[G];				//sm--smoothing function
	double ksqi[G];				//ksqi--1/epsi*(sm/(k*sin(kdx/2)/(kdx/2)))^2

	for(int j=0;j!=G;++j)			//assignment
	{
		k[j]=2*M_PI/(G*dx)*j;
		sm[j]=0;
		ksqi[j]=0;
	}

	for(int j=0;j!=G;++j)			//caculate the smooting function sm and ksqi
	{
		sm[j]=exp(a1*sin(k[j]*dx/2)*sin(k[j]*dx/2)-a2*tan(k[j]*dx/2)*tan(k[j]*dx/2)*tan(k[j]*dx/2)*tan(k[j]*dx/2));
		if(k[j]!=0)
			ksqi[j]=sm[j]*sm[j]/(epsi*k[j]*k[j]*(sin(k[j]*dx/2)/(k[j]*dx/2))*(sin(k[j]*dx/2)/(k[j]*dx/2)));
		else if(k[j]==0)
			ksqi[j]=1.0/(epsi*1);
	}

	Complex *rhok=new Complex[G];		//caculate the rhok[G] via rhoj[G].
	for(int j=0;j!=G;++j)
		rhok[j]=rhoj[t][j];
	CArray data(rhok,G);
	FFT fft;				//fft subroutine
	fft.fft(data);
	for(int j=0;j!=G;++j){
		rhok[j]=data[j];
	}

	Complex *phik=new Complex[G];		//caculate the phik[G] via rhok[G]
	for(int j=0;j!=G;++j){
		phik[j]=rhok[j]*ksqi[j];
		data[j]=phik[j];
	}

	fft.ifft(data);				//caculate the phij[G] via phik[G]
	double *phij=new double[G];
	for(int j=0;j!=G;++j){
		phij[j]=real(data[j]);
	}

	for(int j=0;j!=G;++j)			//caculate the ej[t][128] via phij[t][128].
	{
		if(j==0)
			ej[t][j]=(phij[G-1]-phij[j+1])/(2*dx);
		else if(j==G-1)
			ej[t][j]=(phij[j-1]-phij[0])/(2*dx);
		else
		    ej[t][j]=(phij[j-1]-phij[j+1])/(2*dx);
	}
	
	delete []rhok;				//after the use of the pointer, release the memory on the stack.
	delete []phik;
	delete []phij;
}


