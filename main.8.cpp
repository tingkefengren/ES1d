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

double uniform_dist(double imin,double imax){
	int temp;
	while ((temp=rand()) == RAND_MAX){
	    ;
	 }
	return ((double) temp / RAND_MAX*(imax-imin)+imin);
}

//The main
int main()
{
//Define const variable
/*T--the totle time N--the totle number of the particles G--the totle number of the grids*/
	const int T=500;
	const int N=1000;
	const int G=128;

	double uniform_dist(double imin,double imax);
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

	double q,m,l;
	double dx,dt;
	double epsi,w_p;

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

        /*xi--position of particle
          vi--velocity of particle
          ei--electrostatic of particle
          ese--energy of particle
          p--momenten of particle
          rhoji--charge density of particle
          rhoj--charge density of field
          ej--electrostatic of field
          Xj--position of grid*/

	//Assignment   initial begining.
	q=1.0;
	m=1.0;
	l=0.01;
	epsi=1.0;
	dx=l/G;
	w_p=pow((N/l)*q*q/(epsi*m),0.5);
	dt=0.1*2*M_PI/w_p;

	for(int j=0;j!=G;++j)
		Xj[j]=j*dx;
	for(int t=0;t!=T;++t)
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

	srand((unsigned)time(NULL));
	for(int i=0;i!=N;++i){
		xi[0][i]=uniform_dist(0,1);
		vi[0][i]=2*pow(10,1);
	}

	//Setrho subroutine
	SETRHO rho;
	rho.setrho(xi,Xj,rhoj,q,m,dx,0);

	//Fields subroutine--fields for grid 
	FIELD field;
	field.field(Xj,rhoj,rhoji,ej,epsi,dx,0);

	ofstream E_J("ej.txt");
	for(int j=0;j!=G;++j)
		E_J<<ej[0][j]<<" ";
	E_J<<endl;

        //Field for particle
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

	//Setv subroutine
	ACCEL accel;
	accel.accel(vi,ei,m,q,dt,dx,0);

	//Write the data to the txt
	ofstream X_I("xi.txt");

	/*ES11<<"xi[0][5000]"<<endl;*/
	for(int i=0;i!=N;++i)
	    X_I<<xi[0][i]<<" ";
        X_I<<endl;

	//History subroutine
	HISTORY history;
	history.histry(vi,xi,ese,p,m,0);

        //initial ending

	//The main loop
	for(int t=1;t!=T;++t)
	{   
		//Accel subroutine--accel the particle and change the data of vi[t][5000]
		accel.accel(vi,ei,m,q,dt,dx,t);

		//Move subroutine--change the data of the xi[t][5000]
		MOVE move;
		move.move(vi,xi,dt,dx,t);

		//Write the data of xi[t][5000]and vi[t][5000] at the time t to the txt.

		/*ES11<<"xi[t][5000]"<<t*dt<<endl;*/
		for(int i=0;i!=N;++i)
		    X_I<<xi[t][i]<<" ";
                X_I<<endl;

		//Setrho subroutine--recaculate the electic density of the grid.
		rho.setrho(xi,Xj,rhoj,q,m,dx,t);

		//History subroutine--write the data which is change with the time.
		history.histry(vi,xi,ese,p,m,t);


		//Field subroutine--caculate the electric field intensity via the electric density of the grid.
		field.field(Xj,rhoj,rhoji,ej,epsi,dx,t);

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

	X_I.close();
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
