#include<stdlib.h>
#include<fstream>
#include<iostream>
#define _USE_MATH_DEFINES
#include<math.h>
#include<time.h>
#include "HISTORY.h"
#include "SETRHO.h"
#include "FFT.h"
#include "FIELD.h"
#include "ACCEL.h"
#include "MOVE.h"
using namespace std;

extern const int T=500;		//T--the totle time steps
extern const int N=1000000;	//N--the totle number of the particles
extern const int G=128;		//G--the totle number of the grids

//The uniform distribution
double uniform_dist(double imin,double imax){
	int temp;
	while ((temp=rand()) == RAND_MAX){
	    ;
	 }
	return ((double) temp / RAND_MAX*(imax-imin)+imin);
}

//The maxiwell distribution
double maxwell_dist(double ava,double sig){
	return sig*sqrt(-2.*log(uniform_dist(0,1)))*cos(M_PI*2.*uniform_dist(0,1))+ava;
}

//The main
int main(int argc,char *argv[])
{
	double uniform_dist(double imin,double imax);
	double maxwell_dist(double ava,double sig);
	//Define the pointer to pointer
	double (**xi);			//xi--position of particle
	double (**vi);			//vi--velocity of particle
	double (**ei);                  //ei--electrostatic of particle
	double (**p);                   //p--momenten of particle
	double (**ese);                 //ese--energy of particle
	double *Xj;                     //Xj--position of grid
	double (**rhoj);                //rhoji--charge density of particle
	double (**rhoji);               //rhoj--charge density of field
	double (**ej);                  //ej--electrostatic of field

	int n_0;			//n_0--the number of particle a superparticle certains
	double q,m,epsi;		//q--charge of electron,m--mess of electron,epsi--permittivity
	double l;			//l--length of the system
	double dx,dt;			//dx--length of the grid,dt--time interval
	double x_1,k_mode;		//x_1--amplitude of perturbation,k_mode--mode of perturbation
	double w_p,v_0,v_T;		//w_p--plasma frequency,v_0--drift velocity,v_T--thermal velocity

	//The memory allocation of the pointer to pointer on the stock
	xi=new double*[T];
	vi=new double*[T];
	ei=new double*[T];
	ese=new double*[T];
	p=new double*[T];
	rhoj=new double*[T];
	rhoji=new double*[T];
	ej=new double*[T];
	Xj=new double[G];

	for(int i=0;i!=T;++i)
	{
		xi[i]=new double[N];
		vi[i]=new double[N];
		ei[i]=new double[N];
		ese[i]=new double[N];
		p[i]=new double[N];
		rhoj[i]=new double[G];
		rhoji[i]=new double[G];
		ej[i]=new double[G];
	}

	for(int t=0;t!=T;++t)		//initial the pointer to pointer.
	{
		for(int i=0;i!=N;++i)
		{
			xi[t][i]=0;
			vi[t][i]=0;
			ei[t][i]=0.0;
			ese[t][i]=0.0;
			p[t][i]=0.0;
		}
		for(int j=0;j!=G;++j)
		{
			rhoj[t][j]=0.0;
			rhoji[t][j]=0.0;
			ej[t][j]=0.0;
		}
	}

	////////********initial begin********////////
	n_0=1000;
	q=1.6e-19*n_0;
	m=9.10938215e-31*n_0;
	epsi=8.854187817e-12;
	l=0.01;
	dx=l/G;
	w_p=pow((N/l)*q*q/(epsi*m),0.5);
	dt=0.1*2*M_PI/w_p;
	v_0=0;
	v_T=500;
	x_1=0.1*dx*1000;
	k_mode=0;

	for(int j=0;j!=G;++j)				//initial of grid
		Xj[j]=j*dx;

	srand((unsigned)time(NULL));
	for(int i=0;i!=N;++i){
		xi[0][i]=uniform_dist(0,l);		//initial of xi
		xi[0][i]+=x_1*cos(2*M_PI*k_mode*xi[0][i]/l);	//plus the initial perturbation

		if(xi[0][i]>=l)				//boundary condition
			xi[0][i]-=l;
		else if(xi[0][i]<0)
			xi[0][i]+=l;

		if(v_T==0)				//initial of vi
			vi[0][i]=v_0;
		else
			vi[0][i]=maxwell_dist(0,v_T);
	}

	SETRHO rho;					//Setrho subroutine
	rho.setrho(xi,Xj,rhoj,q,m,dx,0);

	FIELD field;					//Fields subroutine
	field.field(Xj,rhoj,rhoji,ej,epsi,dx,0,argv[1],argv[2]);

	ofstream E_J("ej.txt");
	for(int j=0;j!=G;++j)
		E_J<<ej[0][j]<<" ";
	E_J<<endl;

	for(int i=0;i!=N;++i)				//Field for particle
	{
		for(int j=0;j!=G;++j)
		{
			if(j==G-1)
				if((xi[0][i]>=Xj[j])&&(xi[0][i]<Xj[j]+dx))
					ei[0][i]=((Xj[j]+dx-xi[0][i])/dx)*ej[0][j]+((xi[0][i]-Xj[j])/dx)*ej[0][0];
			else
				if((xi[0][i]>=Xj[j])&&(xi[0][i]<Xj[j+1]))
					ei[0][i]=((Xj[j+1]-xi[0][i])/dx)*ej[0][j]+((xi[0][i]-Xj[j])/dx)*ej[0][j+1];
		}
	}

	ACCEL accel;					//Setv subroutine
	accel.accel(vi,ei,m,q,dt,dx,0);

	ofstream V_I("vi.txt");
	for(int i=0;i!=N;++i)				//Write the data to the txt
		V_I<<vi[0][i]<<" ";
	V_I<<endl;

	ofstream X_I("xi.txt");
	for(int i=0;i!=N;++i)
		X_I<<xi[0][i]<<" ";
        X_I<<endl;

	HISTORY history;				//History subroutine
	history.history(vi,xi,ese,p,m,0);
        ////////********initial end********////////

	////////********The main loop begin********////////
	for(int t=1;t!=T;++t)
	{   
		//Accel subroutine--accel the particle and change the data of vi[t][N]
		accel.accel(vi,ei,m,q,dt,dx,t);

		for(int i=0;i!=N;++i)
			V_I<<vi[t][i]<<" ";
		V_I<<endl;

		//Move subroutine--change the data of the xi[t][N]
		MOVE move;
		move.move(vi,xi,dt,dx,t,l);

		for(int i=0;i!=N;++i)
		    X_I<<xi[t][i]<<" ";
                X_I<<endl;

		//Setrho subroutine--recaculate the electic density of the grid.
		rho.setrho(xi,Xj,rhoj,q,m,dx,t);

		//History subroutine--write the data which is change with the time.
		history.history(vi,xi,ese,p,m,t);

		//Field subroutine--caculate the electric field intensity via the electric density of the grid.
		field.field(Xj,rhoj,rhoji,ej,epsi,dx,t,argv[1],argv[2]);

		for(int j=0;j!=G;++j)
			E_J<<ej[t][j]<<" ";
		E_J<<endl;

		//Caculate the electric field intensity of each particle.
		for(int i=0;i!=N;++i)
		{
			for(int j=0;j!=G;++j)
			{
				if(j==G-1)
					if((xi[0][i]>=Xj[j])&&(xi[0][i]<Xj[j]+dx))
						ei[0][i]=((Xj[j]+dx-xi[0][i])/dx)*ej[0][j]+((xi[0][i]-Xj[j])/dx)*ej[0][0];
				else
					if((xi[0][i]>=Xj[j])&&(xi[0][i]<Xj[j+1]))
						ei[0][i]=((Xj[j+1]-xi[0][i])/dx)*ej[0][j]+((xi[0][i]-Xj[j])/dx)*ej[0][j+1];
			}
		}
	}
	////////********the main loop end********////////

	X_I.close();					//close the file
	V_I.close();
	E_J.close();
	
	//After the use of the pointer, release the memory of the stock.
	for(int i=0;i!=T;++i)
	{
		delete []xi[i];
		delete []vi[i];
		delete []ei[i];
		delete []ese[i];
		delete []p[i];
		delete []rhoj[i];
		delete []rhoji[i];
		delete []ej[i];
	}
	delete []Xj;
	delete []xi;
	delete []vi;
	delete []ei;
	delete []ese;
	delete []p;
	delete []rhoj;
	delete []rhoji;
	delete []ej;

	return 0;
}
