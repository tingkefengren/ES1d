%色散关系分析
%将时域电场信号Ej(x,t)转变为频域谱信号Ek(k,w)
%取波数平均——对谱信号取模，乘以相应的波数并积分（求和），再作归一化

clear
G=128;
T=50;
N=1000;
n_0=1000;
q=1.6e-19*n_0;
m=9.10938215e-31*n_0;
epsi=8.854187817e-12;
l=0.01;
w_p=((N/l)*q*q/(epsi*m))^0.5;
dt=0.1*2*pi/w_p;
t=dt*T;
dx=l/G;
ej=load('ej.txt');
ej_0=sum(ej')/G;
%for i=1:1:T
%	for j=1:1:G
%		ej(i,j)=ej(i,j)-ej_0(i);
%	end
%end
t_i=0:dt:(t-dt);
%t_i=t_i';
x_i=0:dx:(l-dx);
figure;
mesh(x_i,t_i,ej);
xlabel('x(m)');
ylabel('t(s)');
zlabel('E(V/m)');
figure;
contour(x_i,t_i,ej);
xlabel('x(m)');
ylabel('t(s)');

for K=-(G/2):1:(G/2-1)
	k(K+G/2+1)=2*pi*K/l;	%波数
end
for W=-(T/2):1:(T/2-1)
	w(W+T/2+1)=2*pi*W/t;
end
ek=fft2(ej);
ek1=abs(ek);			%取模
figure;
mesh(k,w,ek1);
xlabel('k(m-1)');
ylabel('w(s-1)');
zlabel('功率谱');
figure;
contour(k,w,ek1);
xlabel('k(m-1)');
ylabel('w(s-1)');
%for j=1:1:128
%	ek2_1(j)=0;
%	ek2_2(j)=0;
%	ek2_a(j)=0;
%end
%
%for N=-(G/2):1:(G/2-1)
%	k(N+G/2+1)=2*pi*N/l;	%波数
%end
%
%for j=1:1:128
%	for i=1:1:5000
%		ek2_1(j)=ek2_1(j)+ek1(i,j)*k(j);	%.\frac{\int k dw}{\int dw} .%
%		ek2_a(j)=ek2_a(j)+ek1(i,j);
%	end
%	ek2_2(j)=ek2_1(j)/ek2_a(j);
%end
%plot(k,ek2_2,'.');
%figure;
%plot(k,ek2_a);
