#include<iostream>
#include<stdlib.h>
#include "MOVE.h"
using namespace std;

extern const int N;

void MOVE::move(double (**vi),double (**xi),double dt,double dx,int t,double l)
{
	for(int i=0;i!=N;++i)
	{
		if(t!=0)
		{
			
			xi[t][i]=xi[t-1][i]/dx+vi[t][i]*dt/dx;	//Turn the data to the machine data
			xi[t][i]=xi[t][i]*dx;			//Turn the machine data to the real data
			if(abs(int(xi[t][i])/l)>2)		//Consider the condition of the boundary periodic
				cout<<"Error: the velocity is to big"<<endl;
			else if(xi[t][i]>=l && xi[t][i]<=2*l)
				xi[t][i]=xi[t][i]-l;
			else if(xi[t][i]<0.0 && xi[t][i]>=-2*l)
				xi[t][i]=xi[t][i]+l;
		}
	}
}
