load('xi.txt');
load('vi.txt');
figure;
filename='langmuir.gif'
for i =1:300;
    plot(xi(i,1:100:178000),'.')
    F = getframe(gca);
    im = frame2im(F);
    [I,map] = rgb2ind(im,256);
    k=i-0;
    if k==1;
        imwrite(I,map,filename,'gif','Loopcount',inf, 'DelayTime',0.1);%loopcount只是在i==1的时候才有用
    else
        imwrite(I,map,filename,'gif','WriteMode','append', 'DelayTime',0.1);
    end
	%movie(F)
end
figure;
filename='velocity_hist.gif'
for i =1:300;
    hist(vi(i,:),1000);
    axis([-5 5 0 800]);
    F = getframe(gca);
    im = frame2im(F);
    [I,map] = rgb2ind(im,256);
    k=i-0;
    if k==1;
        imwrite(I,map,filename,'gif','Loopcount',inf, 'DelayTime',0.1);%loopcount只是在i==1的时候才有用
    else
        imwrite(I,map,filename,'gif','WriteMode','append', 'DelayTime',0.1);
    end
	%movie(F)
end
