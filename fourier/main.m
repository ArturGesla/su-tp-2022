

% a=@(n)(-1)^n*(exp(pi)-exp(-pi))/pi/(n^2+1);
% b=@(n)-n*a(n);
% f=a(0)/2;

% a=@(n)((-1)^n-1)/(n^2*pi^2);
% b=@(n)0;
% f=1/4;

a=@(n)0;
b=@(n)(-1)^(n+1)/pi/n;
f=0;


x=[-2*pi:0.01:2*pi];
nm=50;
for n=1:nm
%     f=f+a(n)*cos(n*x)+b(n)*sin(n*x);
    f=f+a(n)*cos(pi/2*n*x)+b(n)*sin(pi/2*n*x);
end
close all;
figure("Position",[100,100,297,210])
plot(x,f); hold on;
% plot(x,exp(x));
plot(x,x/4);
ylim([min(f) max(f)])
xlim([min(x) max(x)]); grid on;
%%
% exportgraphics(gcf,"plot"+num2str(nm)+".png","Resolution",300);
exportgraphics(gcf,"plot5-"+num2str(nm)+".png","Resolution",300);