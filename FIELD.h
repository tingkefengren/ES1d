#ifndef CLASSFIELD_H
#define CLASSFIELD_H

class FIELD
{
public:
	void field(double *Xj,double *rhoj,double *ej,double epsi,double dx,int t,double a1,double a2);
	void get_ei(double *Xj,double *xi,double *ej,double *ei,double dx);
};

#endif
