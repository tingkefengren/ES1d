#include<iostream>
#include "INPUT.h"
using namespace std;

extern "C"
{
	#include "lua.h"
	#include "lauxlib.h"
	#include "lualib.h"
}

extern int T;
extern int G;
extern int N;

void INPUT::lua_to_cxx(Input_p *input_p)
{
	lua_State *L=luaL_newstate();		//creat a lua state
	if(L==NULL){
		//return 0;
		;
	}

	int input=luaL_loadfile(L,"../config/input.lua");	//load the lua file
	if(input){
		cout<<"load file error"<<endl;
	}

	input=lua_pcall(L,0,0,0);		//run the lua file
	if(input){
		cout<<"pcall error"<<endl;
		//return 0;
	}

	lua_getglobal(L,"T");
	T=lua_tonumber(L,-1);
	lua_getglobal(L,"G");
	G=lua_tonumber(L,-1);
	lua_getglobal(L,"N");
	N=lua_tonumber(L,-1);
	lua_getglobal(L,"n");			//read the lua variable
	input_p->n=lua_tonumber(L,-1);
	lua_getglobal(L,"q");
	input_p->q=lua_tonumber(L,-1);
	lua_getglobal(L,"m");
	input_p->m=lua_tonumber(L,-1);
	lua_getglobal(L,"epsi");
	input_p->epsi=lua_tonumber(L,-1);
	lua_getglobal(L,"v_0");
	input_p->v_0=lua_tonumber(L,-1);
	lua_getglobal(L,"v_T");
	input_p->v_T=lua_tonumber(L,-1);
	lua_getglobal(L,"k_mode");
	input_p->k_mode=lua_tonumber(L,-1);
	lua_getglobal(L,"a1");
	input_p->a1=lua_tonumber(L,-1);
	lua_getglobal(L,"a2");
	input_p->a2=lua_tonumber(L,-1);
	lua_getglobal(L,"step_save");
	input_p->step_save=lua_tonumber(L,-1);
	lua_getglobal(L,"normalize");
	input_p->normalize=lua_tonumber(L,-1);
	lua_getglobal(L,"reactive_f_s");
	input_p->reactive_f_s=lua_tonumber(L,-1);
	lua_getglobal(L,"perturb");
	input_p->perturb=lua_tonumber(L,-1);
	lua_getglobal(L,"local_perturb");
	input_p->local_perturb=lua_tonumber(L,-1);
	lua_getglobal(L,"b_position_perturb");
	input_p->b_position_perturb=lua_tonumber(L,-1);
	lua_getglobal(L,"e_position_perturb");
	input_p->e_position_perturb=lua_tonumber(L,-1);
	lua_getglobal(L,"time_to_save_data");
	input_p->time_to_save_data=lua_tonumber(L,-1);
	lua_getglobal(L,"lambda_r");
	input_p->lambda_r=lua_tonumber(L,-1);
	lua_getglobal(L,"w_p_r");
	input_p->w_p_r=lua_tonumber(L,-1);
	lua_getglobal(L,"dt");
	input_p->dt=lua_tonumber(L,-1);
	lua_getglobal(L,"dx");
	input_p->dx=lua_tonumber(L,-1);
	lua_getglobal(L,"l");
	input_p->l=lua_tonumber(L,-1);
	lua_getglobal(L,"x_mode");
	input_p->x_mode=lua_tonumber(L,-1);
	lua_getglobal(L,"n_0");
	input_p->n_0=lua_tonumber(L,-1);
	lua_getglobal(L,"ion_rho");
	input_p->ion_rho=lua_tonumber(L,-1);
	lua_close(L);				//close the lua state
}
