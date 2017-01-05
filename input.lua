T=300			--the total time step
n=1.78e19		--the density
n_0=1e12		--n_0--the number of particle a superparticle certains
q=1.6e-19		--q--charge of electron
m=9.10938215e-31	--m--mess of electron
epsi=8.854187817e-12	--epsi--permittivity
l=0.01			--the length of the system
v_0=0			--the drift velocity
v_T=1.3e6			--the themal velocity
x_mode=0		--the number of times of Perturbation amplitude to dx
k_mode=7		--the wavelength
a1=0			--compensation factor
a2=0			--smoothing factor
step_save=1		--time step to save
normalize=1		--whether to normalize
wp=(n*q^2/(epsi*m))^0.5
lambda=v_T/wp
------------------------------------------
-------------normalization----------------

n_1=n*lambda
q_1=1
m_1=1
epsi_1=epsi*m*lambda*wp^2/q^2
l_1=l/lambda
v_T_1=1

