#include<iostream>
#include "SETRHO.h"
using namespace std;

#define N 1000
#define G 128

void SETRHO::setrho(double (**xi),double *Xj,double (**rhoj),double q,double m,double dx,int t)
{
	//Weighing the charge of the particle to the grid.
	for(int j=0;j!=G;++j)
	{
		for(int i=0;i!=N;++i)
		{
			if(j!=G-1){
				if((xi[t][i]>=Xj[j])&&(xi[t][i]<Xj[j+1])){
					rhoj[t][j]=rhoj[t][j]+q*(Xj[j+1]-xi[t][i])/dx;
					rhoj[t][j+1]=rhoj[t][j+1]+q*(xi[t][i]-Xj[j])/dx;
				}
			}
			else if(j==G-1){
				if((xi[t][i]>=Xj[j])&&(xi[t][i]<Xj[j]+dx)){
					rhoj[t][j]=rhoj[t][j]+q*(Xj[j+1]-xi[t][i])/dx;
					rhoj[t][0]=rhoj[t][0]+q*(xi[t][i]-Xj[j])/dx;
				}
			}
		}
	}
}
