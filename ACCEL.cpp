#include<iostream>
#include "ACCEL.h"
using namespace std;

extern int N;

void ACCEL::accel(double *vi,double *ei,double m,double q,double dt,double dx,int t)
{
	for(int i=0;i!=N;++i)
	{
		
		if(t==0)		//accelerate in the inverse direction about 1/2dt.
			vi[i]=vi[i]-0.5*q*ei[i]*dt/m;
		else			//accelerate in the direct direction about dt.
			vi[i]=vi[i]+q*ei[i]*dt/m;
	}
}
