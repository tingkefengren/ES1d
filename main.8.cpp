#include<stdlib.h>
#include<fstream>
#include<assert.h>
#include<iostream>
#define _USE_MATH_DEFINES
#include<math.h>
#include "HISTORY.h"
#include "SETRHO.h"
#include "FFT.h"
#include "FIELD.h"
#include "ACCEL.h"
#include "MOVE.h"
using namespace std;

//Define variable
//T--the totle time N--the totle number of the particles G--the totle number of the grids
#define T 50
#define N 60
#define G 128

//The main
int main()
{
	//Define the pointer to pointer
	double (**xi);
	double (**vi);
	double (**ei);
	double (**p);
	double (**ese);
	double *Xj;
	double (**rhoj);
	double (**rhoji);
	double (**ej);
	double (**eji);
	double q,m,c;
	double a=3,b=8;
	double dl,dx,dt;
	double epsi;

	//The memory allocation of the pointer to pointer on the stock
	xi=new double*[T];
	vi=new double*[T];
	ei=new double*[T];
	ese=new double*[T];
	p=new double*[T];
	rhoj=new double*[T];
	rhoji=new double*[T];
	ej=new double*[T];
	eji=new double*[T];
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
		eji[i]=new double[G];
	}


	//Assignment
	q=1.0;
	m=1.0;
	epsi=1.0;
	dl=1/127.0;
	dx=1/127.0;
	c=pow(a,b);
	dt=dx/c;
	for(int j=0;j!=G;++j)
		Xj[j]=j*dl;
	for(int a=0;a!=T;++a)
	{
		for(int i=0;i!=N;++i)
		{
			//The particle number density is linear, and the gradient is 1/(5000*5000).
			xi[a][i]=(i+1)/5001.0-(2500.0*i+2500.0-0.5*(1+i*i+2*i))/(5000.0*5000.0);
			vi[a][i]=-0.0003;
			ei[a][i]=0.0;
			ese[a][i]=0.0;
			p[a][i]=0.0;
		}
		for(int j=0;j!=G;++j)
		{
			rhoj[a][j]=0.0;
			rhoji[a][j]=0.0;
			ej[a][j]=0.0;
			eji[a][j]=0.0;
		}
	}

	//Setrho subroutine
	SETRHO rho;
	rho.setrho(xi,Xj,rhoj,q,m,dx,0);

	//Fields subroutine 
	FIELD field;
	field.field(Xj,rhoj,rhoji,ej,eji,epsi,dx,0);

	for(int i=0;i!=N;++i)
	{
		for(int j=0;j!=(G-1);++j)
		{
			if((xi[0][i]>=Xj[j])&&(xi[0][i]<=Xj[j+1]))
			ei[0][i]=((Xj[j+1]-xi[0][i])/dx)*ej[0][j]+((xi[0][i]-Xj[j])/dx)*ej[0][j+1];
		}
	}

	//Setv subroutine
	ACCEL accel;
	accel.accel(vi,ei,m,q,dt,dx,0);

	//Write the data to the txt
	fstream ES11("es1Result8x.txt", fstream::out);
	//fstream ES12("es1Result8v.txt", fstream::out);
	//fstream ES13("es1Result8p.txt", fstream::out);
	//fstream ES14("es1Result8e.txt", fstream::out);
        assert(ES11.is_open());
	//assert(ES12.is_open());
	//assert(ES13.is_open());
	//assert(ES14.is_open());

	//Write the data of xi[0][5000] and vi[0][5000]
	//ES12<<"vi[0][5000]"<<endl;
	//for(int i=0;i!=N;++i)
	//    ES12<<vi[0][i]<<" ";
        //ES12<<endl;

	//ES11<<"xi[0][5000]"<<endl;
	for(int i=0;i!=N;++i)
	    ES11<<xi[0][i]<<" ";
        ES11<<endl;

	//History subroutine
	HISTORY history;
	history.histry(vi,xi,ese,p,m,0);

	//Write the data of ese[0][5000] and p[0][5000]
	//ES14<<"ese[0][5000]"<<endl;
	//for (int i=0;i!=N;++i)
	//	ES14<<ese[0][i]<<endl;

	//ES13<<"p[0][5000]"<<endl;
	//for(int i=0;i!=N;++i)
	//	ES13<<p[0][i]<<endl;

	//The main loop
	for(int t=1;t!=T;++t)
	{   
		//Accel subroutine--accel the particle and change the data of vi[t][5000]
		accel.accel(vi,ei,m,q,dt,dx,t);

		//Move subroutine--change the data of the xi[t][5000]
		MOVE move;
		move.move(vi,xi,dt,dx,t);

		//Write the data of xi[t][5000]and vi[t][5000] at the time t to the txt.
		//ES12<<"vi[t][5000]"<<t*dt<<endl;

		//for(int i=0;i!=N;++i)
		//	ES12<<vi[t][i]<<" ";
                //ES12<<endl;

		//ES11<<"xi[t][5000]"<<t*dt<<endl;
		for(int i=0;i!=N;++i)
		    ES11<<xi[t][i]<<" ";
                ES11<<endl;

		//Setrho subroutine--recaculate the electic density of the grid.
		rho.setrho(xi,Xj,rhoj,q,m,dx,t);

		//History subroutine--write the data which is change with the time.
		history.histry(vi,xi,ese,p,m,t);
        //ES14<<"ese[t][5000]"<<t*dt<<endl;

		//for(int i=0;i!=N;++i)
		//	ES14<<ese[t][i]<<endl;
		
		//ES13<<"p[t][5000]"<<t*dt<<endl;

		//for(int i=0;i!=N;++i)
		//	ES13<<p[t][i]<<endl;

		//Field subroutine--caculate the electric field intensity via the electric density of the grid.
		field.field(Xj,rhoj,rhoji,ej,eji,epsi,dx,t);

		//Caculate the electric field intensity of each particle.
		for(int i=0;i!=N;++i)
		{
			for(int j=0;j!=(G-1);++j)
			{
				if((xi[t][i]>=Xj[j])&&(xi[t][i]<=Xj[j+1]))
				ei[t][i]=((Xj[j+1]-xi[t][i])/dx)*ej[t][j]+((xi[t][i]-Xj[j])/dx)*ej[t][j+1];
			}
		}
	}

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
		delete []eji[i];
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
	delete []eji;

	return 0;
}
