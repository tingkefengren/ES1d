#include<iostream>
#include "ACCEL.h"
using namespace std;

#define N 60

void ACCEL::accel(double (**vi),double (**ei),double m,double q,double dt,double dx,int t)
{
	for(int i=0;i!=N;++i)
	{
		//Turn the data to the mechine data and accelerate in the inverse direction about 1/2dt.
		if(t==0)
			vi[t][i]=vi[t][i]/dx-0.5*q*ei[t][i]*dt/(m*dx);
		//Turn the data to the mechine data and accelerate in the direct direction about dt.
		else
			vi[t][i]=vi[t-1][i]/dx+q*ei[t-1][i]*dt/dx;

		//Turn the mechine data to the real data.
		vi[t][i]=vi[t][i]*dx;
	}
}
