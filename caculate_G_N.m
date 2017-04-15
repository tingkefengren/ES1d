clear;
l=0.01;
n=1.78e16;
n0=1e12;
k_mode=4;
v_T=1.3e6;
me=9.1e-31;
epsi=8.85e-12;
q=1.6e-19;
w_p=sqrt(n*q^2/(me*epsi));
lambda=v_T/w_p;
omega=sqrt(w_p^2+6*pi^2*k_mode^2*v_T^2/l^2);
%dx=0.1*v_T/omega;
dx=lambda;
dt=0.1*2*pi/omega;
%dt=0.1*2*pi/w_p;
G=2^(ceil(log2(l/dx)));
N=n*l/n0;
dx=l/G;
rho=n*q;
e1=rho/epsi*dx;
a=e1*q/me;
dV=a*dt;
dX=dV*dt;
k_max=N/2;
k_min=1;
omega_max=sqrt(w_p^2+6*pi^2*k_max^2*v_T^2/l^2);
omega_min=sqrt(w_p^2+6*pi^2*k_min^2*v_T^2/l^2);
lambda_min=l/k_max;

l_1=l/lambda;
n_1=n*lambda;
v_T_1=1;
me_1=1;
q_1=1;
epsi_1=epsi*me*lambda*w_p^2/q^2;
w_p_1=sqrt(n_1*q_1^2/(me_1*epsi_1));
lambda_1=v_T_1/w_p_1;
omega_1=sqrt(w_p_1^2+6*pi^2*k_mode^2*v_T_1^2/l_1^2);
dx_1=0.1*v_T_1/omega_1;
dt_1=0.1*2*pi/omega_1;
G_1=2^(ceil(log2(l_1/dx_1)));
N_1=n_1*l_1/n0;
dx_1=l_1/G_1;
rho_1=n_1*q_1;
e1_1=rho_1/epsi_1*dx_1;
a_1=e1_1*q_1/me_1;
dV_1=a_1*dx_1;
dX_1=dV_1*dt_1;
