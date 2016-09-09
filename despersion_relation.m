%色散关系分析
%将时域电场信号Ej(x,t)转变为频域谱信号Ek(k,w)
%取波数平均——对谱信号取模，乘以相应的波数并积分（求和），再作归一化

clear
l=0.01;
G=128;
ej=load('ej.txt');
ek=fft2(ej);
ek1=abs(ek);			%取模
for j=1:1:128
	ek2_1(j)=0;
	ek2_2(j)=0;
	ek2_a(j)=0;
end

for N=-(G/2):1:(G/2-1)
	k(N+G/2+1)=2*pi*N/l;	%波数
end

for j=1:1:128
	for i=1:1:50
		ek2_1(j)=ek2_1(j)+ek1(i,j)*k(j);	%.\frac{\int k dw}{\int dw} .%
		ek2_a(j)=ek2_a(j)+ek1(i,j);
	end
	ek2_2(j)=ek2_1(j)/ek2_a(j);
end
plot(k,ek2_2,'.');
figure;
plot(k,ek2_a,'.');
