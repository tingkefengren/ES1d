#include<iostream>
#include<fstream>
#include<vector>
#include<stdlib.h>
#include "INIT.h"
#include "INPUT.h"
#define _USE_MATH_DEFINES
#include<math.h>
using namespace std;

extern int G;
extern int N;

void INIT::initial(vector<double> &Xj,vector<double> &xi,vector<double> &vi,Input_p *input_p)
{
	for(int j=0;j!=G;++j)				//initial of grid
		Xj[j]=j*input_p->dx;
	
	ofstream X_00("../data_analysis/xi00.txt");
	ofstream X_01("../data_analysis/xi01.txt");
	ofstream V_00("../data_analysis/vi00.txt");
	srand((unsigned)time(NULL));

	int ptc_grid=N/G;
	double temp_grid_f=0;						//a temporary grid former
	double temp_grid_l=temp_grid_f+input_p->dx;						//a temporary grid latter

	//for(int i=0;i!=N-N%G;++i){
	//	xi[i]=uniform_dist(temp_grid_f,temp_grid_l);		//initial of xi
	//	if((i+1)%ptc_grid==0){
	//		temp_grid_f+=dx_r;
	//		temp_grid_l+=dx_r;
	//	}
	//}
	//for(int i=0;i!=N%G;++i)
	//	xi[i+N-N%G]=uniform_dist(0,l_r);

	double v_T_beam=0.3*input_p->v_T;
	for(int i=0;i!=N;++i)
		xi[i]=uniform_dist(0,input_p->l);		//initial of xi

	if(input_p->local_perturb==1){
		for(int i=0;i!=N;++i){
			if(xi[i]>input_p->b_position_perturb && xi[i]<input_p->e_position_perturb)
				xi[i]+=input_p->x_mode*cos(2*M_PI*input_p->k_mode*xi[i]/input_p->l);	//plus the initial perturbation
			if(xi[i]>=input_p->l)				//boundary condition
				xi[i]-=input_p->l;
			else if(xi[i]<0)
				xi[i]+=input_p->l;
			X_01<<xi[i]<<" ";

			if(input_p->v_T==0)				//initial of vi
				vi[i]=input_p->v_0;
			else{
				if(i%20==0)
					vi[i]=maxwell_dist(input_p->v_0,v_T_beam);
				else
					vi[i]=maxwell_dist(0,input_p->v_T);
			}
			V_00<<vi[i]<<" ";
		}
	}
	else if(input_p->local_perturb==0){
		for(int i=0;i!=N;++i){
			X_00<<xi[i]<<" ";
			xi[i]+=input_p->x_mode*cos(2*M_PI*input_p->k_mode*xi[i]/input_p->l);	//plus the initial perturbation

			if(xi[i]>=input_p->l)				//boundary condition
				xi[i]-=input_p->l;
			else if(xi[i]<0)
				xi[i]+=input_p->l;
			X_01<<xi[i]<<" ";

			if(input_p->v_T==0)				//initial of vi
				vi[i]=input_p->v_0;
			else{
				if(i%20==0)
					vi[i]=maxwell_dist(input_p->v_0,v_T_beam);
				else
					vi[i]=maxwell_dist(0,input_p->v_T);
			}
			V_00<<vi[i]<<" ";

		}
	}
	X_00.close();
	X_01.close();
	V_00.close();
}

double INIT::uniform_dist(double imin,double imax)
{
	int temp;
	while ((temp=rand()) == RAND_MAX){
	    ;
	 }
	return ((double) temp / RAND_MAX*(imax-imin)+imin);
}

double INIT:: maxwell_dist(double mean,double sigma)
{
	double ymin=mean-4.*sigma;
	double ymax=mean+4.*sigma;
	double Pymax=1./sqrt(2.*M_PI)/sigma;
	double y=ymin+(ymax-ymin)*double(random())/double(RAND_MAX);
	double Py=exp(-(y-mean)*(y-mean)/2./sigma/sigma)/sqrt(2.*M_PI)/sigma;
	double x=Pymax*double(random())/double(RAND_MAX);
	if(x>Py)
		return maxwell_dist(mean,sigma);
	else
		return y;
}

//The maxiwell distribution
//double INIT::maxwell_dist(double ava,double sig){
	//return sig*sqrt(-2.*log(uniform_dist(0,1)))*cos(M_PI*2.*uniform_dist(0,1))+ava;
	
//}


