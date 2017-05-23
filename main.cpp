#include<stdlib.h>
#include<iomanip>
#include<fstream>
#include<vector>
#include<iostream>
#define _USE_MATH_DEFINES
#include<math.h>
#include<mpi.h>
#include<time.h>
#include "HISTORY.h"
#include "SETRHO.h"
#include "FIELD.h"
#include "ACCEL.h"
#include "MOVE.h"
#include "INIT.h"
#include "INPUT.h"
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


//The main
int main(int argc,char* argv[])
{
	Input_p input_p;
	INPUT INput_p;
	INput_p.lua_to_cxx(&input_p);
	////****get data from lua file end****////

	vector<double> xi(N);			//xi--position of particle
	vector<double> vi(N);			//vi--velocity of particle
	vector<double> ei(N);                  //ei--electrostatic of particle
	vector<double> p(N);                   //p--momenten of particle
	vector<double> ese(N);                 //ese--energy of particle
	vector<double> Xj(G);                     //Xj--position of grid
	vector<double> rhoj(G);               //rhoj--charge density of field
	vector<double> ej(G);                  //ej--electrostatic of field
	
	int thread_size,my_rank,group_size;
	MPI::Init(argc,argv);
	my_rank=MPI::COMM_WORLD.Get_rank();
	thread_size=MPI::COMM_WORLD.Get_size();
	group_size=N/thread_size;

	INIT init;
	if(my_rank==0){
		init.initial(Xj,xi,vi,&input_p);
		cout<<"The Particle Number N: "<<N<<setw(20)<<" "<<"The Total Time Step T: "<<T<<endl;
		cout<<"The Time Interval dt: "<<input_p.dt<<"s"<<setw(15)<<" "<<"The Grid Interval dx: "<<input_p.dx<<"m"<<endl;
	}

	////////********initial begin********////////
	if(input_p.normalize==1){
		if(my_rank==0){
			for(int i=0;i!=N;++i){
				xi[i]=xi[i]/input_p.lambda_r;
				vi[i]=vi[i]/(input_p.lambda_r*input_p.w_p_r);
			}
			for(int j=0;j!=G;++j){
				Xj[j]=Xj[j]/input_p.lambda_r;
			}
		}
	}
	

	vector<double> rhoj_mid(G);

	MPI::COMM_WORLD.Bcast(&xi.front(),xi.size(),MPI::DOUBLE,0);
	MPI::COMM_WORLD.Bcast(&Xj.front(),Xj.size(),MPI::DOUBLE,0);

	SETRHO rho;					//Setrho subroutine
	rho.setrho(my_rank,group_size,xi,Xj,rhoj,input_p.q,input_p.n_0,input_p.dx);		//get recent rhoj[1][j]. rhoj[0][j]=0
	MPI::COMM_WORLD.Barrier();

	if(my_rank!=0){
		MPI::COMM_WORLD.Send(&rhoj.front(),rhoj.size(),MPI::DOUBLE,0,0);
	}
	else{
		for(int j=0;j!=G;++j)
			rhoj_mid[j]=rhoj[j];
		for(int i=1;i!=thread_size;++i){
			MPI::COMM_WORLD.Recv(&rhoj.front(),rhoj.size(),MPI::DOUBLE,i,0);
			if(i!=thread_size-1){
				for(int j=0;j!=G;++j)
					rhoj_mid[j]+=rhoj[j];
			}
			else if(i==thread_size-1){
				for(int j=0;j!=G;++j)
					rhoj[j]=rhoj_mid[j]+rhoj[j];
			}
		}
	}

	if(my_rank==0){
		for(int j=0;j!=G;++j){
			rhoj[j]=(rhoj[j]+input_p.ion_rho)/input_p.dx;
		}
	ofstream RHO("../data_analysis/rho00.txt");
	for(int j=0;j!=G;++j)
		RHO<<rhoj[j]<<" ";
	RHO.close();
	}

	ofstream E_J("../data_analysis/ej.txt");
	ofstream V_I("../data_analysis/vi.txt");
	ofstream X_I("../data_analysis/xi.txt");
	FIELD field;					//Fields subroutine
	int t=0;
	if(my_rank==0){
		field.field(Xj,rhoj,ej,input_p.epsi,input_p.dx,t,input_p.a1,input_p.a2);	//get recent ej[1][j]. ej[0][j]=0

		if(t>=input_p.time_to_save_data){
			for(int j=0;j!=G;++j)
				E_J<<ej[j]<<" ";
			E_J<<endl;
		}
	}

	MPI::COMM_WORLD.Bcast(&ej.front(),ej.size(),MPI::DOUBLE,0);

	field.get_ei(my_rank,group_size,Xj,xi,ej,ei,input_p.dx);			//Field for particle
	MPI::COMM_WORLD.Barrier();
	vector<double> ei_mid(N);

	if(my_rank!=0){
		MPI::COMM_WORLD.Send(&ei.front(),ei.size(),MPI::DOUBLE,0,0);
	}
	else{
		for(int i=0;i!=group_size;++i)
			ei_mid[i]=ei[i];
		for(int i=1;i!=thread_size;++i){
			MPI::COMM_WORLD.Recv(&ei.front(),ei.size(),MPI::DOUBLE,i,0);
			for(int k=i*group_size;k!=(i+1)*group_size;++k)
				ei_mid[k]=ei[k];
		}
		for(int i=0;i!=N;++i)
			ei[i]=ei_mid[i];
	}

	ACCEL accel;					//Setv subroutine
	HISTORY history;				//History subroutine
	//MPI::COMM_WORLD.Bcast(ei,N,MPI::DOUBLE,0);
	//accel.accel(my_rank,group_size,vi,ei,m,q,dt,dx,0);
	//MPI::COMM_WORLD.Barrier();
	//double* vi_mid=new double[N];


	if(my_rank==0){
		accel.accel(vi,ei,input_p.m,input_p.q,input_p.dt,input_p.dx,t);

		if(t>=input_p.time_to_save_data){
			for(int i=0;i!=N;++i)				//Write the data to the txt
				V_I<<vi[i]<<" ";
			V_I<<endl;

			for(int i=0;i!=N;++i)
				X_I<<xi[i]<<" ";
    	    X_I<<endl;
		}

		history.history(vi,xi,ese,p,input_p.m,t);
        ////////********initial end********////////
	}
	////////********The main loop begin********////////
	for(int t=1;t!=T;++t)
	{   
		if(my_rank==0){
			//Accel subroutine--accel the particle and change the data of vi[t][N]
			accel.accel(vi,ei,input_p.m,input_p.q,input_p.dt,input_p.dx,t);

			//Move subroutine--change the data of the xi[t][N]
			MOVE move;
			move.move(vi,xi,input_p.dt,input_p.dx,t,input_p.l);

			if(t>=input_p.time_to_save_data){
				if(input_p.step_save!=0){
					if(t%input_p.step_save==0){
						for(int i=0;i!=N;++i)
							V_I<<vi[i]<<" ";
						V_I<<endl;

						for(int i=0;i!=N;++i)
						    X_I<<xi[i]<<" ";
        			        	X_I<<endl;
					}
				}
			}
		}
		MPI::COMM_WORLD.Bcast(&xi.front(),xi.size(),MPI::DOUBLE,0);
		//Setrho subroutine--recaculate the electic density of the grid.
		for(int j=0;j!=G;++j){
			rhoj[j]=0;
			rhoj_mid[j]=0;
		}
		//clear the data of rhoj[0][j] in previous step

		rho.setrho(my_rank,group_size,xi,Xj,rhoj,input_p.q,input_p.n_0,input_p.dx);
		MPI::COMM_WORLD.Barrier();

		if(my_rank!=0){
			MPI::COMM_WORLD.Send(&rhoj.front(),rhoj.size(),MPI::DOUBLE,0,0);
		}
		else{
			for(int j=0;j!=G;++j)
				rhoj_mid[j]=rhoj[j];
			for(int i=1;i!=thread_size;++i){
				MPI::COMM_WORLD.Recv(&rhoj.front(),rhoj.size(),MPI::DOUBLE,i,0);
				if(i!=thread_size-1){
					for(int j=0;j!=G;++j)
						rhoj_mid[j]+=rhoj[j];
				}
				else if(i==thread_size-1){
					for(int j=0;j!=G;++j)
						rhoj[j]=rhoj_mid[j]+rhoj[j];
				}
			}
		}

		if(my_rank==0){
			for(int j=0;j!=G;++j){
				rhoj[j]=(rhoj[j]+input_p.ion_rho)/input_p.dx;
			}

			//History subroutine--write the data which is change with the time.
			history.history(vi,xi,ese,p,input_p.m,t);

			//Field subroutine--caculate the electric field intensity via the electric density of the grid.
			field.field(Xj,rhoj,ej,input_p.epsi,input_p.dx,t,input_p.a1,input_p.a2);

			if(t>=input_p.time_to_save_data){
				if(input_p.step_save!=0){
					if(t%input_p.step_save==0){
						for(int j=0;j!=G;++j)
							E_J<<ej[j]<<" ";
						E_J<<endl;
					}
				}
			}

			//Caculate the electric field intensity of each particle.
		}

		MPI::COMM_WORLD.Bcast(&ej.front(),ej.size(),MPI::DOUBLE,0);

		for(int i=0;i!=N;i++){
			ei_mid[i]=0.0;
			ei[i]=0.0;
		}
			
		field.get_ei(my_rank,group_size,Xj,xi,ej,ei,input_p.dx);
		MPI::COMM_WORLD.Barrier();

		if(my_rank!=0){
			MPI::COMM_WORLD.Send(&ei.front(),ei.size(),MPI::DOUBLE,0,0);
		}
		else{
			for(int i=0;i!=group_size;++i)
				ei_mid[i]=ei[i];
			for(int i=1;i!=thread_size;++i){
				MPI::COMM_WORLD.Recv(&ei.front(),ei.size(),MPI::DOUBLE,i,0);
				for(int k=i*group_size;k!=(i+1)*group_size;++k)
					ei_mid[k]=ei[k];
				}
			for(int i=0;i!=N;++i)
				ei[i]=ei_mid[i];
		}
	}

	////////********the main loop end********////////

	MPI::Finalize();
	X_I.close();					//close the file
	V_I.close();
	E_J.close();
	
	return 0;
}
