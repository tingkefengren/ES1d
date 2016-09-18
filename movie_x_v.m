%getframe
clear
xi=load('xi.txt');
vi=load('vi.txt');
count=1;
for i=1:1:50

%    contour(1:33, 1:33, real(phi(:,:,i)),20);

	plot(xi(i,:),vi(i,:),'.');

	axesValue = axis ;

	axis(axesValue) ;

	M(count)=getframe(gcf);

	count = count + 1;

end

% making movie
myObj = VideoWriter('trace_x_v.avi');%初始化一个avi文件

myObj.FrameRate = 5;

myObj.Quality = 100;

open(myObj);
for i=1:count-1


	writeVideo(myObj,M(i));
end
close(myObj);


