#include<iostream>
#include<stdlib.h>
#include "MOVE.h"
using namespace std;

extern const int N;

void MOVE::move(double *vi,double *xi,double dt,double dx,int t,double l)
{
	for(int i=0;i!=N;++i)
	{
		if(t!=0)
		{
			xi[i]=xi[i]/dx+vi[i]*dt/dx;	//Turn the data to the machine data
			xi[i]=xi[i]*dx;			//Turn the machine data to the real data
			if(abs(int(xi[i])/l)>2)		//Consider the condition of the boundary periodic
				cout<<"Error: the velocity is to big"<<endl;
			else if(xi[i]>=l && xi[i]<=2*l)
				xi[i]=xi[i]-l;
			else if(xi[i]<0.0 && xi[i]>=-2*l)
				xi[i]=xi[i]+l;
		}
	}
}
