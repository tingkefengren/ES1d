%色散关系分析
%将时域电场信号Ej(x,t)转变为频域谱信号Ek(k,w)
%取波数平均——对谱信号取模，乘以相应的波数并积分（求和），再作归一化

clear;
T=300;
l_r=0.01;
n_r=1.78e19;
n_0=1e12;
q_r=1.6e-19;
m_r=9.10938215e-31;
epsi_r=8.854187817e-12;
v_T_r=1e6;
k_mode=7;
step_save=10;

w_p_r=(n_r*q_r^2/(epsi_r*m_r))^0.5;
lambda_r=v_T_r/w_p_r;
omega_r=(w_p_r^2+1.5*(2*pi*k_mode/l_r)^2*v_T_r^2)^0.5;
dt_r=0.1/w_p_r;
t_r=dt_r*T;
dx_r=0.1*v_T_r/omega_r;
G_r=2^(ceil(log2(l_r/dx_r)));
dx_r=l_r/G_r;

q=1;
m=1;
v_T=1;
l=l_r/lambda_r;
n=n_r*lambda_r;
epsi=epsi_r*m_r*lambda_r*w_p_r^2/q_r^2;
w_p=(n*q^2/(epsi*m))^0.5;
omega=(w_p^2+1.5*(2*pi*k_mode/l)^2*v_T^2)^0.5;
dt=0.1/w_p;
t=dt*T;
dx=0.1*v_T/omega;
G=2^(ceil(log2(l/dx)));
dx=l/G;

t_s=dt:dt:T*dt;

ej=load('ej.txt');
ej_e=0.5*sum((ej.^2*dx*epsi)');
t_i=0:dt:(t-dt);
x_i=0:dx:(l-dx);

figure;%('visible','off')
plot(t_i,ej_e);
xlabel('t(s)');
ylabel('electric field energy(J)');
%saveas(gcf,'elctric_energy.eps','epsc');

vi=load('vi.txt');
motion_e=0.5*m*n_0.*sum((vi.^2)');
figure;%('visible','off')
plot(t_i,motion_e);
xlabel('t(s)');
ylabel('motion energy(J)');
%saveas(gcf,'motion_energy.eps','epsc');

figure;%('visible','off')
total_e=motion_e+ej_e;
plot(t_i,total_e);
xlabel('t(s)');
ylabel('total energy(J)');
%saveas(gcf,'total_energy.eps','epsc');

figure;%('visible','off');
mesh(x_i,t_i,ej);
xlabel('x(m)');
ylabel('t(s)');
zlabel('E(V/m)');
%saveas(gcf,'ej_mesh.eps','epsc');

figure;%('visible','off');
contourf(x_i,t_i,ej);
xlabel('x(m)');
ylabel('t(s)');
%saveas(gcf,'ej_contour.eps','epsc');

for K=0:1:G-1
	k(K+1)=2*pi*K/l;	%波数
end
for W=0:1:T-1
	w(W+1)=2*pi*W/t;	%频率
end
ek=fft2(ej);
ek1=abs(ek);			%取模
ek2=ek1(1:T/2,1:G/2);

figure;%('visible','off');
mesh(k(1:G/2),w(1:T/2),ek2);
xlabel('k(m-1)');
ylabel('w(s-1)');
zlabel('功率谱');
%saveas(gcf,'ek_mesh.eps','epsc');

figure;%('visible','off');
contourf(k(1:G/2),w(1:T/2),ek2);
xlabel('k(m-1)');
ylabel('w(s-1)');
colormap(cool);

figure;%('visible','off');
[minx,omega_n]=min(abs(w-omega));
plot(k(1:G/2),ek2(omega_n,1:G/2));
xlabel('k');
ylabel('功率谱');
title('朗缪尔波功率谱');
%saveas(gcf,'ek_mesh.eps','epsc');
