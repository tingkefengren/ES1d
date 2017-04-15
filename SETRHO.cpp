#include<iostream>
#include "SETRHO.h"
using namespace std;

extern int N;
extern int G;

void SETRHO::setrho(double *xi,double *Xj,double *rhoj,double q,double n_0,double dx,int t,double ion_rho)
{
	//Weighing the charge of the particle to the grid.
	int j=0;
	for(int i=0;i!=N;++i)
	{
		j=(int)(xi[i]/dx);
		if(j!=G-1){
				rhoj[j]=rhoj[j]+q*n_0*(Xj[j+1]-xi[i])/dx;
				rhoj[j+1]=rhoj[j+1]+q*n_0*(xi[i]-Xj[j])/dx;
		}
		else if(j==G-1){
				rhoj[j]=rhoj[j]+q*n_0*(Xj[j]+dx-xi[i])/dx;
				rhoj[0]=rhoj[0]+q*n_0*(xi[i]-Xj[j])/dx;
		}
	}
	for(int j=0;j!=G;++j){
		rhoj[j]=(rhoj[j]+ion_rho)/dx;
	}
}
