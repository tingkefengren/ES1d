#ifndef CLASSMOVE_H
#define CLASSMOVE_H
#include<vector>
using namespace std;

class MOVE
{
public:
	void move(vector<double> &vi,vector<double> &xi,double &dt,double &dx,int &t,double &l);
	//void move(double *vi,double *xi,double dt,double dx,int t,double l);
};

#endif
