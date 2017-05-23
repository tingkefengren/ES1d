#include<iostream>
#include<stdlib.h>
#include<vector>
#include "MOVE.h"
using namespace std;

extern int N;

void MOVE::move(vector<double> &vi,vector<double> &xi,double &dt,double &dx,int &t,double &l)
//void MOVE::move(double *vi,double *xi,double dt,double dx,int t,double l)
{
	for(int i=0;i!=N;++i)
	{
		if(t!=0)
		{
			xi[i]=xi[i]+vi[i]*dt;	//Turn the data to the machine data
			if(abs(int(xi[i])/l)>2)		//Consider the condition of the boundary periodic
				cout<<"Error: the velocity is to big"<<endl;
			else if(xi[i]>=l && xi[i]<=2*l)
				xi[i]=xi[i]-l;
			else if(xi[i]<0.0 && xi[i]>=-2*l)
				xi[i]=xi[i]+l;
		}
	}
}
