T=200			--the total time step
G=2048
N=178000
n=1.78e13		--the density
q=-1.6e-19		--q--charge of electron
m=9.10938215e-31	--m--mess of electron
epsi=8.854187817e-12	--epsi--permittivity
v_0=0			--the drift velocity
v_T=1e6			--the themal velocity
k_mode=3		--the wavelength
a1=0			--compensation factor
a2=0			--smoothing factor
step_save=1		--time step to save
normalize=0		--whether to normalize
reactive_f_s=0	--choose for reactive problem for dt, 1 means the fast mode, 0 means the slow mode
local_perturb=0		--whether to add the perturbation locally
b_position_perturb=0.2	--the begining position of perturbation
e_position_perturb=0.8	--the ending posion of perturbation
wp=(n*q^2/(epsi*m))^0.5
lambda=v_T/wp
------------------------------------------
-------------normalization----------------

n_1=n*lambda
q_1=1
m_1=1
epsi_1=epsi*m*lambda*wp^2/q^2

