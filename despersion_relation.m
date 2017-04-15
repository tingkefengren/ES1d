%色散关系分析
%将时域电场信号Ej(x,t)转变为频域谱信号Ek(k,w)
%取波数平均——对谱信号取模，乘以相应的波数并积分（求和），再作归一化

clear;
T=200;
G=2048;
N=178000;
n_r=1.78e13;
q_r=1.6e-19;
m_r=9.10938215e-31;
epsi_r=8.854187817e-12;
v_T_r=1e6;
v_0=0;
k_mode=100;
step_save=1;
normalize=0;

w_p_r=(n_r*q_r^2/(epsi_r*m_r))^0.5;
if v_T_r==0
	lambda_r=0.02*v_0/w_p_r;
else
	lambda_r=v_T_r/w_p_r;
end
dx_r=lambda_r;
l_r=G*dx_r;
x_mode=l_r/(40*k_mode);
n_0=n_r*l_r/N;
omega_r=(w_p_r^2+1.5*(2*pi*k_mode/l_r)^2*v_T_r^2)^0.5;
dt_r=0.05*2*pi/w_p_r;
t_r=dt_r*T;
dx_r=lambda_r;

if normalize==1
	q=1;
	m=1;
	v_T=1;
	l=l_r/lambda_r;
	n=n_r*lambda_r;
	epsi=epsi_r*m_r*lambda_r*w_p_r^2/q_r^2;
	w_p=(n*q^2/(epsi*m))^0.5;
	omega=(w_p^2+1.5*(2*pi*k_mode/l)^2*v_T^2)^0.5;
	dt=0.05*2*pi/w_p;
	t=dt*T;
	dx=v_T/w_p;
	G=2^(ceil(log2(l/dx)));
	dx=l/G;
elseif normalize==0
	q=q_r;
	m=m_r;
	v_T=v_T_r;
	l=l_r;
	n=n_r;
	epsi=epsi_r;
	w_p=w_p_r;
	omega=omega_r;
	dt=dt_r;
	t=dt*T;
	dx=dx_r;
end

ej=load('ej.txt');
ej_e=0.5*sum((ej.^2*dx*epsi)');
t_i=0:dt:(t-dt);
x_i=0:dx:(l-dx);

figure;%('visible','off')
plot(t_i,ej_e);
xlabel({'t(s)'},'Interpreter','latex');
ylabel({'electric field energy(J)'},'Interpreter','latex');
title({'Electric Field Energy Change with Time'},'Interpreter','latex');
%saveas(gcf,'elctric_energy.eps','epsc');

%vi=load('vi.txt');
%motion_e=0.5*m*n_0.*sum((vi.^2)');
%figure;%('visible','off')
%plot(t_i,motion_e);
%xlabel({'t(s)'},'Interpreter','latex');
%ylabel({'particle motion energy(J)'},'Interpreter','latex');
%title({'Particle Motion Energy Change with Time'},'Interpreter','latex');
%%%saveas(gcf,'motion_energy.eps','epsc');
%%
%figure;%('visible','off')
%total_e=motion_e+ej_e;
%plot(t_i,total_e);
%ylabel('total energy(J)');
%xlabel({'t(s)'},'Interpreter','latex');
%ylabel({'total energy(J)'},'Interpreter','latex');
%title({'Total Energy Change with Time'},'Interpreter','latex');
%%saveas(gcf,'total_energy.eps','epsc');
%
figure;%('visible','off');
mesh(x_i,t_i,ej);
%surf(x_i,t_i,ej);
xlabel({'x(m)'},'Interpreter','latex');
ylabel({'t(s)'},'Interpreter','latex');
zlabel({'E(V/m)'},'Interpreter','latex');
title({'Electirc Field Change with Space and Time'},'Interpreter','latex');
%%saveas(gcf,'ej_mesh.eps','epsc');

figure;%('visible','off');
contourf(x_i,t_i,ej);
xlabel({'x(m)'},'Interpreter','latex');
ylabel({'t(s)'},'Interpreter','latex');
title({'Electirc Field Change with Space and Time'},'Interpreter','latex');
%saveas(gcf,'ej_contour.eps','epsc');

f=2*pi/(T*dt)*(-T/2:T/2-1);
k=2*pi/(G*dx)*(-G/2:G/2-1);

ej=ej-repmat(mean(ej,2),1,G);
ek=fft2(ej);
ek1=fftshift(abs(ek));			%取模

figure;%('visible','off');
%mesh(k1,w,ek2);
mesh(k,f,ek1);
%surf(k2,w,ek2);
xlabel({'$k(m^{-1})$'},'Interpreter','latex');
ylabel({'$f (s^{-1})$'},'Interpreter','latex');
zlabel({'power spectrum'},'Interpreter','latex');
title({'The Power Spectrum of Electirc Field'},'Interpreter','latex');
%%saveas(gcf,'ek_mesh.eps','epsc');
%
figure;%('visible','off');
%contourf(k1,w,ek2);
contourf(k,f,ek1);
xlabel({'$k(m^{-1})$'},'Interpreter','latex');
ylabel({'$f (s^{-1})$'},'Interpreter','latex');
title({'The Power Spectrum of Electirc Field'},'Interpreter','latex');
colormap(cool);

figure;%('visible','off');
[minx,omega_n]=min(abs(f-omega));
%plot(k1,ek2(omega_n,1:G/2+1));
plot(k,ek1(omega_n,:));
xlabel({'$k(m^{-1})$'},'Interpreter','latex');
ylabel({'power spectrum'},'Interpreter','latex');
title({'$The\ Power\ Spectrum\ of\ Electirc\ Field\ at\ frenquence\ \omega$'},'Interpreter','latex');
%saveas(gcf,'ek_mesh.eps','epsc');
figure;%('visible','off');
[minx,k_n]=min(abs(k-2*pi*k_mode/l));
plot(f,ek1(:,k_n-1));
%%title('朗缪尔波功率谱');
xlabel({'$\omega \ (s^{-1})$'},'Interpreter','latex');
ylabel({'power spectrum'},'Interpreter','latex');
title({'$The\ Power\ Spectrum\ of\ Electirc\ Field\ at\ mode\ k\_ mode$'},'Interpreter','latex');

