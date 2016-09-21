#include<iostream>
#include "ACCEL.h"
using namespace std;

extern const int N;

void ACCEL::accel(double (**vi),double (**ei),double m,double q,double dt,double dx,int t)
{
	for(int i=0;i!=N;++i)
	{
		
		if(t==0)		//accelerate in the inverse direction about 1/2dt.
			vi[t][i]=vi[t][i]/dx-0.5*q*ei[t][i]*dt/(m*dx);
		else			//accelerate in the direct direction about dt.
			vi[t][i]=vi[t-1][i]/dx+q*ei[t-1][i]*dt/(m*dx);

		vi[t][i]=vi[t][i]*dx;	//Turn the mechine data to the real data.
	}
}
