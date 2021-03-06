#include<iostream>
#include<vector>
#include<mpi.h>
#include "SETRHO.h"
using namespace std;

extern int N;
extern int G;

void SETRHO::setrho(int &my_rank,int &group_size,vector<double> &xi,vector<double> &Xj,vector<double> &rhoj,double &q,double &n_0,double &dx)
//void SETRHO::setrho(int my_rank,int group_size,double *xi,double *Xj,double *rhoj,double q,double n_0,double dx)
{


	//Weighing the charge of the particle to the grid.
	int j=0;
	for(int i=my_rank*group_size;i!=(my_rank+1)*group_size;++i)
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

}
