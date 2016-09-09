#include<iostream>
#include "MOVE.h"
using namespace std;

#define N 1000

void MOVE::move(double (**vi),double (**xi),double dt,double dx,int t)
{
	for(int i=0;i!=N;++i)
	{
		if(t==0)
		{
			//Turn the data to the mechine data.
			xi[t][i]=xi[t][i]/dx;
			//Turn the mechine data to the real data.
			xi[t][i]=xi[t][i]*dx;
		}
		else
		{
			//Turn the data to the machine data.
			xi[t][i]=xi[t-1][i]/dx+vi[t][i]*dt/dx;
			//Turn the machine data to the real data.
			xi[t][i]=xi[t][i]*dx;
			//Consider the condition of the boundary periodic.
			if(xi[t][i]>1.0)
				xi[t][i]=xi[t][i]-int(xi[t][i]);
			else if(xi[t][i]<0.0)
				xi[t][i]=int(xi[t][i])-xi[t][i];
		}
	}
}
