#ifndef CLASSFIELD_H
#define CLASSFIELD_H
#include<vector>
using namespace std;

class FIELD
{
public:
	void field_fft(vector<double> &Xj,vector<double> &rhoj,vector<double> &ej,double &epsi,double &dx,int &t,double &a1,double &a2);
	void field_hockey(vector<double> &Xj,vector<double> &rhoj,vector<double> &ej,double &epsi,double &dx,int &t,double &a1,double &a2);
	void get_ei(int &my_rank,int &group_size,vector<double> &Xj,vector<double> &xi,vector<double> &ej,vector<double> &ei,double &dx);
	//void field(double *Xj,double *rhoj,double *ej,double epsi,double dx,int t,double a1,double a2);
	//void get_ei(int my_rank,int group_size,double *Xj,double *xi,double *ej,double *ei,double dx);
};

#endif
