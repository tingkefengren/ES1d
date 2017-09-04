#ifndef CLASSINPUT_H
#define CLASSINPUT_H
using namespace std;

typedef struct INPUT_P
{
	double n,q,m,epsi,v_0,v_T,k_mode,a1,a2,a3;
	int step_save,normalize,reactive_f_s,perturb,local_perturb;
	int field_model;
	double b_position_perturb,e_position_perturb,time_to_save_data,lambda_r,w_p_r,dt,dx,l,x_mode,n_0,ion_rho,different_dist_v_instability,linear_dist,slope_linear_dist,b_linear_dist,e_linear_dist; 
}Input_p;

class INPUT
{
public:
	void lua_to_cxx(Input_p *input_p);
};
#endif
