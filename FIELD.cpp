#include<iostream>
#include<fstream>
#include<vector>
#include<complex>
#include<valarray>
#include<stdlib.h>
#include "FIELD.h"
#include "FFT.h"
#define _USE_MATH_DEFINES
#include<math.h>
using namespace std;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

extern int G;
extern int N;

void FIELD::field_fft(vector<double> &Xj,vector<double> &rhoj,vector<double> &ej,double &epsi,double &dx,int &t,double &a1,double &a2)
//void FIELD::field(double *Xj,double *rhoj,double *ej,double epsi,double dx,int t,double a1,double a2)
{
	double k[G];				//k--wave number
	double sm[G];				//sm--smoothing function
	double ksqi[G];				//ksqi--1/epsi*(sm/(k*sin(kdx/2)/(kdx/2)))^2

	for(int j=0;j!=G;++j)			//assignment
	{
		if(j<G/2)
			k[j]=2*M_PI/(G*dx)*j;
		else
			k[j]=2*M_PI/(G*dx)*(j-G);
		sm[j]=0;
		ksqi[j]=0;
	}

	for(int j=0;j!=G;++j)			//caculate the smooting function sm and ksqi
	{
		//sm[j]=exp(a1*sin(k[j]*dx/2)*sin(k[j]*dx/2)-a2*tan(k[j]*dx/2)*tan(k[j]*dx/2)*tan(k[j]*dx/2)*tan(k[j]*dx/2));
		if(k[j]!=0)
			//ksqi[j]=sm[j]*sm[j]/(epsi*k[j]*k[j]*(sin(k[j]*dx/2)/(k[j]*dx/2))*(sin(k[j]*dx/2)/(k[j]*dx/2)));
			ksqi[j]=1.0/(epsi*k[j]*k[j]*(sin(k[j]*dx/2)/(k[j]*dx/2))*(sin(k[j]*dx/2)/(k[j]*dx/2)));
		else if(k[j]==0)
			ksqi[j]=1.0/(epsi*1);
	}

	Complex *rhok=new Complex[G];		//caculate the rhok[G] via rhoj[G].
	for(int j=0;j!=G;++j)
		rhok[j]=rhoj[j];
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
			ej[j]=(phij[G-1]-phij[j+1])/(2*dx);
		else if(j==G-1)
			ej[j]=(phij[j-1]-phij[0])/(2*dx);
		else
		    ej[j]=(phij[j-1]-phij[j+1])/(2*dx);
	}
	
	delete []rhok;				//after the use of the pointer, release the memory on the stack.
	delete []phik;
	delete []phij;
}

void FIELD::field_hockey(vector<double> &Xj,vector<double> &rhoj,vector<double> &ej,double &epsi,double &dx,int &t,double &a1,double &a2)
{
	//double *phij=new double[G];
	vector<double> phij(G);
	for(int j=0;j!=G;j++)
		phij[j]=0.0;
	double sum_rho=0;
	for(int j=1;j!=G;++j)
		sum_rho+=j*dx*dx/epsi*(-rhoj[G-j]);
	phij[1]=-1.0/G*sum_rho;
	for(int j=2;j!=G;j++)
		phij[j]=-rhoj[j]*dx*dx/epsi+2*phij[j-1]-phij[j-2];
	for(int j=0;j!=G;++j)			//caculate the ej[t][128] via phij[t][128].
	{
		if(j==0)
			ej[j]=(phij[G-1]-phij[j+1])/(2*dx);
		else if(j==G-1)
			ej[j]=(phij[j-1]-phij[0])/(2*dx);
		else
		    ej[j]=(phij[j-1]-phij[j+1])/(2*dx);
	}
	//delete []phij;
}

void FIELD::get_ei(int &my_rank,int &group_size,vector<double> &Xj,vector<double> &xi,vector<double> &ej,vector<double> &ei,double &dx)
//void FIELD::get_ei(int my_rank,int group_size,double *Xj,double *xi,double *ej,double *ei,double dx)
{
	int j=0;
	for(int i=my_rank*group_size;i!=(my_rank+1)*group_size;++i){
		j=(int)(xi[i]/dx);
		if(j==G-1)
			ei[i]=((Xj[j]+dx-xi[i])/dx)*ej[j]+((xi[i]-Xj[j])/dx)*ej[0];
		else
			ei[i]=((Xj[j+1]-xi[i])/dx)*ej[j]+((xi[i]-Xj[j])/dx)*ej[j+1];
	}
}

