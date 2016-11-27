#include<iostream>
#include "HISTORY.h"
using namespace std;

extern const int N;

void HISTORY::history(double *vi,double *xi,double *ese,double *p,double m,int t)
{
    for(int i=0;i!=N;++i)
	{
		//Caculate the enegetic and the momentum of the paticle.
		ese[i]=0.5*m*vi[i]*vi[i];
		p[i]=m*vi[i];
	}
}
