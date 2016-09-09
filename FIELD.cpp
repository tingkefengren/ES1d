#include<iostream>
#include<fstream>
#include "FIELD.h"
#include "FFT.h"
#define _USE_MATH_DEFINES
#include<math.h>
using namespace std;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

#define G 128

void FIELD::field(double *Xj,double (**rhoj),double (**rhoji),double (**ej),double epsi,double dx,int t)
{
	//Define the variable.
	FFT fft;
	double ng=128.0;
	double kdx2[G];
	double sm[G];
	double ksqi[G];
	Complex *rhok=new Complex[G];
	for(int j=0;j!=G;++j)
		rhok[j]=rhoj[t][j];
	CArray data(rhok,G);
		
	//Assignment
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
	//Caculate the rhok[50][128] via rhoj[50][128].
	//FFT_1D() subroutine
	fft.fft(data);
	for(int j=0;j!=G;++j)
		rhok[j]=data[j];

	//Caculate the phik[t][128] and phiki[t][128] via rhok[t][128] and rhoki[t][128].
	Complex *phik=new Complex[G];
	

	//IFFT_1D subroutine
	fft.ifft(data);

	double *phij=new double[G];
	for(int j=0;j!=G;++j)
		phij[j]=real(data[j]);

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


