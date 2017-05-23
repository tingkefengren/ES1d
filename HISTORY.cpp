#include<iostream>
#include<vector>
#include "HISTORY.h"
using namespace std;

extern int N;

void HISTORY::history(vector<double> &vi,vector<double> &xi,vector<double> &ese,vector<double> &p,double &m,int &t)
//void HISTORY::history(double *vi,double *xi,double *ese,double *p,double m,int t)
{
    for(int i=0;i!=N;++i)
	{
		//Caculate the enegetic and the momentum of the paticle.
		ese[i]=0.5*m*vi[i]*vi[i];
		p[i]=m*vi[i];
	}
}
