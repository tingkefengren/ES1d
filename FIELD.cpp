#include<iostream>
#include<fstream>
#include "FIELD.h"
#include "FFT.h"
#define _USE_MATH_DEFINES
#include<math.h>
using namespace std;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

extern const int G;

void FIELD::field(double *Xj,double (**rhoj),double (**rhoji),double (**ej),double epsi,double dx,int t)
{
	//Define the variable.
	FFT fft;
	double a1,a2;
	double k[G];
	double sm[G];
	double ksqi[G];
	Complex *rhok=new Complex[G];
	for(int j=0;j!=G;++j)
		rhok[j]=rhoj[t][j];
	CArray data(rhok,G);

	//Assignment
	for(int j=0;j!=G;++j)
	{
		k[j]=(2*M_PI/(G*dx))*(j-G/2);
		sm[j]=0;
		ksqi[j]=0;
	}

	a1=0;
	a2=0;
	//Caculate the smooting factor sm and 1/(K*K).
	for(int j=0;j!=G;++j)
	{
		sm[j]=exp(a1*sin(k[j]*dx/2)*sin(k[j]*dx/2)-a2*tan(k[j]*dx/2)*tan(k[j]*dx/2)*tan(k[j]*dx/2)*tan(k[j]*dx/2));
		if(k[j]!=0)
			ksqi[j]=sm[j]*sm[j]/(epsi*k[j]*k[j]*(sin(k[j]*dx/2)/(k[j]*dx/2))*(sin(k[j]*dx/2)/(k[j]*dx/2)));
		else if(k[j]==0)
			ksqi[j]=1.0/(epsi*1);
	}
	//Caculate the rhok[50][128] via rhoj[50][128].
	//FFT_1D() subroutine
	fft.fft(data);
	for(int j=0;j!=G;++j){
		rhok[j]=data[j];
	}

	//Caculate the phik[t][128] and phiki[t][128] via rhok[t][128] and rhoki[t][128].
	Complex *phik=new Complex[G];

	for(int j=0;j!=G;++j){
		phik[j]=rhok[j]*ksqi[j];
		data[j]=phik[j];
	}

	//IFFT_1D subroutine
	fft.ifft(data);

	double *phij=new double[G];
	for(int j=0;j!=G;++j){
		phij[j]=real(data[j]);
	}

	//Caculate the ej[t][128] via phij[t][128].
	for(int j=0;j!=G;++j)
	{
		if(j==0)
			ej[t][j]=(phij[G-1]-phij[j+1])/(2*dx);
		else if(j==G-1)
			ej[t][j]=(phij[j-1]-phij[0])/(2*dx);
		else
		    ej[t][j]=(phij[j-1]-phij[j+1])/(2*dx);
	}
	
//	cout<<t;
	//After the use of the pointer, release the memory on the stack.
	delete []rhok;
	delete []phik;
	delete []phij;
}


