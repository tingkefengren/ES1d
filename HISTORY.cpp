#include<iostream>
#include "HISTORY.h"
using namespace std;

#define N 60

void HISTORY::histry(double (**vi),double (**xi),double (**ese),double (**p),double m,int t)
{
    for(int i=0;i!=N;++i)
	{
		//Caculate the enegetic and the momentum of the paticle.
		ese[t][i]=0.5*m*vi[t][i]*vi[t][i];
		p[t][i]=m*vi[t][i];
	}
}
