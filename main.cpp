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
//double maxwell_dist(double ava,double sig){
	//return sig*sqrt(-2.*log(uniform_dist(0,1)))*cos(M_PI*2.*uniform_dist(0,1))+ava;
	
//}

double maxwell_dist(double mean,double sigma)
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
//The main
int main()
{
	double uniform_dist(double imin,double imax);
	double maxwell_dist(double mean,double sigma);
	//Define the pointer to pointer
	double *xi;			//xi--position of particle
	double *vi;			//vi--velocity of particle
	double *ei;                  //ei--electrostatic of particle
	double *p;                   //p--momenten of particle
	double *ese;                 //ese--energy of particle
	double *Xj;                     //Xj--position of grid
	double *rhoj;               //rhoj--charge density of field
	double *ej;                  //ej--electrostatic of field
	double ion_rho_r;				//ion_rho_r--the rho contribute by ion at each grid

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
	lua_getglobal(L,"G");
	G=lua_tonumber(L,-1);
	lua_getglobal(L,"N");
	N=lua_tonumber(L,-1);
	lua_getglobal(L,"n");			//read the lua variable
	double n_r=lua_tonumber(L,-1);
	lua_getglobal(L,"q");
	double q_r=lua_tonumber(L,-1);
	lua_getglobal(L,"m");
	double m_r=lua_tonumber(L,-1);
	lua_getglobal(L,"epsi");
	double epsi_r=lua_tonumber(L,-1);
	lua_getglobal(L,"v_0");
	double v_0_r=lua_tonumber(L,-1);
	lua_getglobal(L,"v_T");
	double v_T_r=lua_tonumber(L,-1);
	lua_getglobal(L,"k_mode");
	double k_mode=lua_tonumber(L,-1);
	lua_getglobal(L,"a1");
	double a1=lua_tonumber(L,-1);
	lua_getglobal(L,"a2");
	double a2=lua_tonumber(L,-1);
	lua_getglobal(L,"step_save");
	int step_save=lua_tonumber(L,-1);
	lua_getglobal(L,"normalize");
	int normalize=lua_tonumber(L,-1);
	lua_getglobal(L,"reactive_f_s");
	int reactive_f_s=lua_tonumber(L,-1);
	lua_getglobal(L,"local_perturb");
	int local_perturb=lua_tonumber(L,-1);
	lua_getglobal(L,"b_position_perturb");
	double b_position_perturb=lua_tonumber(L,-1);
	lua_getglobal(L,"e_position_perturb");
	double e_position_perturb=lua_tonumber(L,-1);
	lua_close(L);				//close the lua state
	////****get data from lua file end****////
	clock_t end_lua=clock();
	TIME_EACH<<(end_lua-start_lua)/60000000<<"min"<<(end_lua-start_lua)%60000000/1000000<<"s"<<(end_lua-start_lua)%60000000%1000000/1000<<"ms"<<(end_lua-start_lua)%60000000%1000000%1000<<"mus"<<endl;

	double w_p_r=pow(n_r*q_r*q_r/(epsi_r*m_r),0.5);
	double lambda_r;
	if(v_T_r!=0)
		lambda_r=v_T_r/w_p_r;
	else if(v_T_r==0)
		lambda_r=0.01*v_0_r/w_p_r;
	//omega_r for electron
	double dx_r=lambda_r;
	double l_r=dx_r*G;
	double x_mode=0.5*l_r/(40*k_mode);
	double n_0=n_r*l_r/N;
	double omega_r=pow((w_p_r*w_p_r+1.5*4*M_PI*M_PI*k_mode*k_mode/(l_r*l_r)*v_T_r*v_T_r),0.5);
	double omega_f_r=2*M_PI*k_mode/l_r*v_0_r+w_p_r;
	double omega_s_r=2*M_PI*k_mode/l_r*v_0_r-w_p_r;
	double dt_r;
	dt_r=0.05*2*M_PI/w_p_r;
		//if(reactive_f_s==1)
		//	dt_r=0.04*2*M_PI/omega_f_r;
		//else if(reactive_f_s==0)
		//	dt_r=0.005*2*M_PI/omega_s_r;
	ion_rho_r=-q_r*n_r*dx_r;

	//The memory allocation of the pointer to pointer on the stock
	xi=new double[N];
	vi=new double[N];
	ei=new double[N];
	ese=new double[N];
	p=new double[N];
	rhoj=new double[G];
	ej=new double[G];
	Xj=new double[G];

	for(int i=0;i!=N;++i)
	{
		xi[i]=0;
		vi[i]=0;
		ei[i]=0.0;
		ese[i]=0.0;
		p[i]=0.0;
	}
	for(int j=0;j!=G;++j)
	{
		Xj[j]=0.0;
		rhoj[j]=0.0;
		ej[j]=0.0;
	}
	

	////////********initial begin********////////
	for(int j=0;j!=G;++j)				//initial of grid
		Xj[j]=j*dx_r;
	
	ofstream X_00("xi00.txt");
	ofstream X_01("xi01.txt");
	ofstream V_00("vi00.txt");
	srand((unsigned)time(NULL));

	int ptc_grid=N/G;
	double temp_grid_f=0;						//a temporary grid former
	double temp_grid_l=temp_grid_f+dx_r;						//a temporary grid latter
	for(int i=0;i!=N-N%G;++i){
		xi[i]=uniform_dist(temp_grid_f,temp_grid_l);		//initial of xi
		if((i+1)%ptc_grid==0){
			temp_grid_f+=dx_r;
			temp_grid_l+=dx_r;
		}
	}
	for(int i=0;i!=N%G;++i)
		xi[i+N-N%G]=uniform_dist(0,l_r);

	if(local_perturb==1){
		for(int i=0;i!=N;++i){
			if(xi[i]>b_position_perturb && xi[i]<e_position_perturb)
				xi[i]+=x_mode*cos(2*M_PI*k_mode*xi[i]/l_r);	//plus the initial perturbation
			if(xi[i]>=l_r)				//boundary condition
				xi[i]-=l_r;
			else if(xi[i]<0)
				xi[i]+=l_r;
			X_01<<xi[i]<<" ";

			if(v_T_r==0)				//initial of vi
				vi[i]=v_0_r;
			else
				vi[i]=maxwell_dist(0,v_T_r);
			V_00<<vi[i]<<" ";
		}
	}
	else if(local_perturb==0){
		for(int i=0;i!=N;++i){
			X_00<<xi[i]<<" ";
			xi[i]+=x_mode*cos(2*M_PI*k_mode*xi[i]/l_r);	//plus the initial perturbation

			if(xi[i]>=l_r)				//boundary condition
				xi[i]-=l_r;
			else if(xi[i]<0)
				xi[i]+=l_r;
			X_01<<xi[i]<<" ";

			if(v_T_r==0)				//initial of vi
				vi[i]=v_0_r;
			else
				vi[i]=maxwell_dist(0,v_T_r);
			V_00<<vi[i]<<" ";

		}
	}
	X_00.close();
	X_01.close();
	V_00.close();

	double l,n,v_T,m,q,epsi,w_p,omega,omega_f,omega_s,v_0,dt,dx,x_1,ion_rho;
	if(normalize==1){
		//////******normalization begin******//////
		l=l_r/lambda_r;
		n=n_r*lambda_r;
		v_T=1;
		v_0=v_0_r/v_T_r;
		m=1;
		q=1;
		epsi=epsi_r*m_r*lambda_r*w_p_r*w_p_r/(q_r*q_r);
		w_p=pow(n*q*q/(epsi*m),0.5);
		//omega for eletron
		omega=pow((w_p*w_p+1.5*4*M_PI*M_PI*k_mode*k_mode/(l*l)*v_T*v_T),0.5);
		omega_f=2*M_PI*k_mode/l*v_0+w_p;
		omega_s=2*M_PI*k_mode/l*v_0-w_p;
		//omega for ion
		//omega=v_T*k_mode;
		dt=0.05*2*M_PI/w_p;
			//if(reactive_f_s==1)
			//	dt=0.04*2*M_PI/omega_f;
			//else if(reactive_f_s==0)
			//	dt=0.005*2*M_PI/omega_s;
		dx=v_T/w_p;
		G=pow(2,ceil(log(l/dx)/log(2)));
		dx=l/G;
		ion_rho=-q*n*dx;

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
		v_0=v_0_r;
		m=m_r;
		q=q_r;
		epsi=epsi_r;
		w_p=w_p_r;
		omega=omega_r;
		omega_f=omega_f_r;
		omega_s=omega_s_r;
		dt=dt_r;
		dx=dx_r;
		ion_rho=ion_rho_r;
	}

	TIME_EACH<<"time of set rho is: ";
	clock_t start_setrho=clock();
	SETRHO rho;					//Setrho subroutine
	rho.setrho(xi,Xj,rhoj,q,n_0,dx,0,ion_rho);		//get recent rhoj[1][j]. rhoj[0][j]=0
	ofstream RHO("rho00.txt");
	for(int j=0;j!=G;++j)
		RHO<<rhoj[j]<<" ";
	RHO.close();

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
	cout<<ei[100]<<endl;
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

		rho.setrho(xi,Xj,rhoj,q,n_0,dx,t,ion_rho);

		//History subroutine--write the data which is change with the time.
		history.history(vi,xi,ese,p,m,t);

		//Field subroutine--caculate the electric field intensity via the electric density of the grid.
		field.field(Xj,rhoj,ej,epsi,dx,t,a1,a2);

		if(t%step_save==0){
			for(int j=0;j!=G;++j)
				E_J<<ej[j]<<" ";
			E_J<<endl;
		}

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
