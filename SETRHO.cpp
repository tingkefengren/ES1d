#include<iostream>
#include "SETRHO.h"
using namespace std;

extern const int N;
extern const int G;

void SETRHO::setrho(double *xi,double *Xj,double *rhoj,double q,double m,double dx,int t)
{
	//Weighing the charge of the particle to the grid.
	for(int j=0;j!=G;++j)
	{
		for(int i=0;i!=N;++i)
		{
			if(j!=G-1){
				if((xi[i]>=Xj[j])&&(xi[i]<Xj[j+1])){
					rhoj[j]=rhoj[j]+q*(Xj[j+1]-xi[i])/dx;
					rhoj[j+1]=rhoj[j+1]+q*(xi[i]-Xj[j])/dx;
				}
			}
			else if(j==G-1){
				if((xi[i]>=Xj[j])&&(xi[i]<Xj[j]+dx)){
					rhoj[j]=rhoj[j]+q*(Xj[j]+dx-xi[i])/dx;
					rhoj[0]=rhoj[0]+q*(xi[i]-Xj[j])/dx;
				}
			}
		}
	}
}
