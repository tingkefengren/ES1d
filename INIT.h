#ifndef CLASSINIT_H
#define CLASSINIT_H
#include "INPUT.h"
#include<vector>
using namespace std;

class INIT
{
public:
	void initial(vector<double> &Xj,vector<double> &xi,vector<double> &vi,Input_p *input_p);
	//void initial(vector<double> &Xj,vector<double> &xi,vector<double> &vi,double &dx,double &l,double &x_mode,double &k_mode,int &local_perturb,double &b_position_perturb,double &e_position_perturb,double &v_0,double &v_T);
	//void initial(double *Xj,double *xi,double *vi,double dx_r,double l_r,double x_mode,double k_mode,int local_perturb,double b_position_perturb,double e_position_perturb,double v_0_r,double v_T_r);
	double uniform_dist(double imin,double imax);
	double maxwell_dist(double mean,double sigma);
	//double INIT::maxwell_dist(double ava,double sig)
};

#endif
