#ifndef CLASSSETRHO_H
#define CLASSSETRHO_H
#include<vector>
using namespace std;

class SETRHO
{
public:
	void setrho(int &my_rank,int &group_size,vector<double> &xi,vector<double> &Xj,vector<double> &rhoj,double &q,double &n_0,double &dx);
	//void setrho(int my_rank,int group_size,double *xi,double *Xj,double *rhoj,double q,double n_0,double dx);
};

#endif
