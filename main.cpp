#include<stdlib.h>
#include<fstream>
#include<iostream>
#define _USE_MATH_DEFINES
#include<math.h>
#include<time.h>
#include "HISTORY.h"
#include "SETRHO.h"
#include "FIELD.h"
#include "ACCEL.h"
#include "MOVE.h"
using namespace std;

extern "C"
{
	#include "lua.h"
	#include "lauxlib.h"
	#include "lualib.h"
}

int T;		//T--the totle time steps
int N;		//N--the totle number of the particles
int G;		//G--the totle number of the grids

//The uniform distribution
double uniform_dist(double imin,double imax){
	int temp;
	while ((temp=rand()) == RAND_MAX){
	    ;
	 }
	return ((double) temp / RAND_MAX*(imax-imin)+imin);
}

//The maxiwell distribution
double maxwell_dist(double ava,double sig){
	return sig*sqrt(-2.*log(uniform_dist(0,1)))*cos(M_PI*2.*uniform_dist(0,1))+ava;
}

//The main
int main()
{
	double uniform_dist(double imin,double imax);
	double maxwell_dist(double ava,double sig);
	//Define the pointer to pointer
	double *xi;			//xi--position of particle
	double *vi;			//vi--velocity of particle
	double *ei;                  //ei--electrostatic of particle
	double *p;                   //p--momenten of particle
	double *ese;                 //ese--energy of particle
	double *Xj;                     //Xj--position of grid
	double *rhoj;               //rhoj--charge density of field
	double *ej;                  //ej--electrostatic of field

	ofstream TIME_EACH("time_of_each_part.txt");
	TIME_EACH<<"time of lua is: ";
	clock_t start_lua=clock();
	////****get data from lua file begin****////
	lua_State *L=luaL_newstate();		//creat a lua state
	if(L==NULL){
		return 0;
	}

	int input=luaL_loadfile(L,"input.lua");	//load the lua file
	if(input){
		cout<<"load file error"<<endl;
	}

	input=lua_pcall(L,0,0,0);		//run the lua file
	if(input){
		cout<<"pcall error"<<endl;
		return 0;
	}

	lua_getglobal(L,"T");
	T=lua_tonumber(L,-1);
	lua_getglobal(L,"n");			//read the lua variable
	double n_r=lua_tonumber(L,-1);
	lua_getglobal(L,"q");
	double q_r=lua_tonumber(L,-1);
	lua_getglobal(L,"m");
	double m_r=lua_tonumber(L,-1);
	lua_getglobal(L,"epsi");
	double epsi_r=lua_tonumber(L,-1);
	lua_getglobal(L,"l");
	double l_r=lua_tonumber(L,-1);
	lua_getglobal(L,"v_0");
	double v_0=lua_tonumber(L,-1);
	lua_getglobal(L,"v_T");
	double v_T_r=lua_tonumber(L,-1);
	lua_getglobal(L,"x_mode");
	double x_mode=lua_tonumber(L,-1);
	lua_getglobal(L,"k_mode");
	double k_mode=lua_tonumber(L,-1);
	lua_getglobal(L,"n_0");
	double n_0=lua_tonumber(L,-1);
	lua_getglobal(L,"a1");
	double a1=lua_tonumber(L,-1);
	lua_getglobal(L,"a2");
	double a2=lua_tonumber(L,-1);
	lua_getglobal(L,"step_save");
	int step_save=lua_tonumber(L,-1);
	lua_getglobal(L,"normalize");
	int normalize=lua_tonumber(L,-1);
	lua_close(L);				//close the lua state
	////****get data from lua file end****////
	clock_t end_lua=clock();
	TIME_EACH<<(end_lua-start_lua)/60000000<<"min"<<(end_lua-start_lua)%60000000/1000000<<"s"<<(end_lua-start_lua)%60000000%1000000/1000<<"ms"<<(end_lua-start_lua)%60000000%1000000%1000<<"mus"<<endl;

	N=n_r*l_r/n_0;
	double w_p_r=pow(n_r*q_r*q_r/(epsi_r*m_r),0.5);
	double lambda_r=v_T_r/w_p_r;
	double omega_r=pow((w_p_r*w_p_r+1.5*4*M_PI*M_PI*k_mode*k_mode/(l_r*l_r)*v_T_r*v_T_r),0.5);
	double dt_r=0.1/w_p_r;
	double dx_r=0.1*v_T_r/omega_r;
	int G_r=pow(2,ceil(log(l_r/dx_r)/log(2)));
	dx_r=l_r/G_r;
	double x_1_r=x_mode*dx_r;

	//The memory allocation of the pointer to pointer on the stock
	xi=new double[N];
	vi=new double[N];
	ei=new double[N];
	ese=new double[N];
	p=new double[N];
	rhoj=new double[G_r];
	ej=new double[G_r];
	Xj=new double[G_r];

	for(int i=0;i!=N;++i)
	{
		xi[i]=0;
		vi[i]=0;
		ei[i]=0.0;
		ese[i]=0.0;
		p[i]=0.0;
	}
	for(int j=0;j!=G_r;++j)
	{
		Xj[j]=0.0;
		rhoj[j]=0.0;
		ej[j]=0.0;
	}
	

	////////********initial begin********////////
	for(int j=0;j!=G_r;++j)				//initial of grid
		Xj[j]=j*dx_r;
	
	srand((unsigned)time(NULL));
	for(int i=0;i!=N;++i){
		xi[i]=uniform_dist(0,l_r);		//initial of xi
		xi[i]+=x_1_r*cos(2*M_PI*k_mode*xi[i]/l_r);	//plus the initial perturbation

		if(xi[i]>=l_r)				//boundary condition
			xi[i]-=l_r;
		else if(xi[i]<0)
			xi[i]+=l_r;


		if(v_T_r==0)				//initial of vi
			vi[i]=v_0;
		else
			vi[i]=maxwell_dist(0,v_T_r);

	}

	double l,n,v_T,m,q,epsi,w_p,omega,dt,dx,x_1;
	if(normalize==1){
		//////******normalization begin******//////
		l=l_r/lambda_r;
		n=n_r*lambda_r;
		v_T=1;
		m=1;
		q=1;
		epsi=epsi_r*m_r*lambda_r*w_p_r*w_p_r/(q_r*q_r);
		w_p=pow(n*q*q/(epsi*m),0.5);
		omega=pow((w_p*w_p+1.5*4*M_PI*M_PI*k_mode*k_mode/(l*l)*v_T*v_T),0.5);
		dt=0.1/omega;
		dx=0.1*v_T/omega;
		G=pow(2,ceil(log(l/dx)/log(2)));
		dx=l/G;
		x_1=x_mode*dx;

		for(int i=0;i!=N;++i){
			xi[i]=xi[i]/lambda_r;
			vi[i]=vi[i]/(lambda_r*w_p_r);
		}
		for(int j=0;j!=G;++j)
			Xj[j]=Xj[j]/lambda_r;
		//////******normalization end******//////
	}
	else if(normalize==0){
		l=l_r;
		n=n_r;
		v_T=v_T_r;
		m=m_r;
		q=m_r;
		epsi=epsi_r;
		w_p=w_p_r;
		omega=omega_r;
		dt=dt_r;
		dx=dx_r;
		G=G_r;
		x_1=x_1_r;
	}

	TIME_EACH<<"time of set rho is: ";
	clock_t start_setrho=clock();
	SETRHO rho;					//Setrho subroutine
	rho.setrho(xi,Xj,rhoj,q,n_0,dx,0);		//get recent rhoj[1][j]. rhoj[0][j]=0
	clock_t end_setrho=clock();
	TIME_EACH<<(end_setrho-start_setrho)/60000000<<"min"<<(end_setrho-start_setrho)%60000000/1000000<<"s"<<(end_setrho-start_setrho)%60000000%1000000/1000<<"ms"<<(end_setrho-start_setrho)%60000000%1000000%1000<<"mus"<<endl;

	TIME_EACH<<"time of field is: ";
	clock_t start_field=clock();
	FIELD field;					//Fields subroutine
	field.field(Xj,rhoj,ej,epsi,dx,0,a1,a2);	//get recent ej[1][j]. ej[0][j]=0
	clock_t end_field=clock();
	TIME_EACH<<(end_field-start_field)/60000000<<"min"<<(end_field-start_field)%60000000/1000000<<"s"<<(end_field-start_field)%60000000%1000000/1000<<"ms"<<(end_field-start_field)%60000000%1000000%1000<<"mus"<<endl;

	ofstream E_J("ej.txt");
	for(int j=0;j!=G;++j)
		E_J<<ej[j]<<" ";
	E_J<<endl;

	TIME_EACH<<"time of get_ei is: ";
	clock_t start_get_ei=clock();
	field.get_ei(Xj,xi,ej,ei,dx);			//Field for particle
	clock_t end_get_ei=clock();
	TIME_EACH<<(end_get_ei-start_get_ei)/60000000<<"min"<<(end_get_ei-start_get_ei)%60000000/1000000<<"s"<<(end_get_ei-start_get_ei)%60000000%1000000/1000<<"ms"<<(end_get_ei-start_get_ei)%60000000%1000000%1000<<"mus"<<endl;

	TIME_EACH<<"time of accel is: ";
	clock_t start_accel=clock();
	ACCEL accel;					//Setv subroutine
	accel.accel(vi,ei,m,q,dt,dx,0);
	clock_t end_accel=clock();
	TIME_EACH<<(end_accel-start_accel)/60000000<<"min"<<(end_accel-start_accel)%60000000/1000000<<"s"<<(end_accel-start_accel)%60000000%1000000/1000<<"ms"<<(end_accel-start_accel)%60000000%1000000%1000<<"mus"<<endl;

	ofstream V_I("vi.txt");
	for(int i=0;i!=N;++i)				//Write the data to the txt
		V_I<<vi[i]<<" ";
	V_I<<endl;

	ofstream X_I("xi.txt");
	for(int i=0;i!=N;++i)
		X_I<<xi[i]<<" ";
        X_I<<endl;

	HISTORY history;				//History subroutine
	history.history(vi,xi,ese,p,m,0);
        ////////********initial end********////////

	TIME_EACH<<"time of main loop is: ";
	clock_t start_main_loop=clock();
	////////********The main loop begin********////////
	for(int t=1;t!=T;++t)
	{   
		//Accel subroutine--accel the particle and change the data of vi[t][N]
		accel.accel(vi,ei,m,q,dt,dx,t);

		//Move subroutine--change the data of the xi[t][N]
		MOVE move;
		move.move(vi,xi,dt,dx,t,l);

		if(t%step_save==0){
			for(int i=0;i!=N;++i)
				V_I<<vi[i]<<" ";
			V_I<<endl;

			for(int i=0;i!=N;++i)
			    X_I<<xi[i]<<" ";
                	X_I<<endl;
		}
		//Setrho subroutine--recaculate the electic density of the grid.
		for(int j=0;j!=G;++j)
			rhoj[j]=0;
		//clear the data of rhoj[0][j] in previous step

		rho.setrho(xi,Xj,rhoj,q,n_0,dx,t);

		//History subroutine--write the data which is change with the time.
		history.history(vi,xi,ese,p,m,t);

		//Field subroutine--caculate the electric field intensity via the electric density of the grid.
		field.field(Xj,rhoj,ej,epsi,dx,t,a1,a2);

		for(int j=0;j!=G;++j)
			E_J<<ej[j]<<" ";
		E_J<<endl;

		//Caculate the electric field intensity of each particle.
		field.get_ei(Xj,xi,ej,ei,dx);
	}
	////////********the main loop end********////////
	clock_t end_main_loop=clock();
	TIME_EACH<<(end_main_loop-start_main_loop)/60000000<<"min"<<(end_main_loop-start_main_loop)%60000000/1000000<<"s"<<(end_main_loop-start_main_loop)%60000000%1000000/1000<<"ms"<<(end_main_loop-start_main_loop)%60000000%1000000%1000<<"mus"<<endl;
	TIME_EACH<<"G: "<<G<<endl;
	TIME_EACH<<"T: "<<T<<endl;

	X_I.close();					//close the file
	V_I.close();
	E_J.close();
	
	//After the use of the pointer, release the memory of the stock.
	delete []Xj;
	delete []xi;
	delete []vi;
	delete []ei;
	delete []ese;
	delete []p;
	delete []rhoj;
	delete []ej;

	return 0;
}
