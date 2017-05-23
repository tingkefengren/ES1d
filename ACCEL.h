#ifndef CLASSACCEL_H
#define CLASSACCEL_H
#include<vector>
using namespace std;

class ACCEL
{
public:
	//void accel(int my_rank,int group_size,double *vi,double *ei,double m,double q,double dt,double dx,int t);
	void accel(vector<double> &vi,vector<double> &ei,double &m,double &q,double &dt,double &dx,int &t);
	//void accel(double *vi,double *ei,double m,double q,double dt,double dx,int t);
};

#endif
