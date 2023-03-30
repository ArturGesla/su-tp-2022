l=2;
a=@(n)0;
% b=@(n)(2/n/pi)*(1-cos(n*(pi)));
f=0;
omega=@(n)(2*n+1)/2*pi/l;
b=@(n)cos(omega(n)*l/2)*(4/l/omega(n)-2/omega(n))+sin(omega(n)*l/2)*4/l/omega(n)^2 ...
+cos(omega(n)*l)*(4/l/omega(n))+sin(omega(n)*l)*(-2)/l/omega(n)^2;

% b=@(n)1/omega(n)*(2/l/omega(n)*sin(omega(n)*l/2)-cos(omega(n)*l/2))+ ...
%     4/l/omega(n)*(cos(omega(n)*l/2)-cos(omega(n)*l)) ...
%     -1/omega(n)*cos(l*omega(n)/2) ...
%     -2/l/omega(n)^2*(sin(omega(n)*l)-sin(omega(n)*l/2));
% x=[-2*pi:0.01:2*pi];
x=[-3*l:0.025:3*l];
nm=50;
for n=0:nm
%     f=f+a(n)*cos(n*x)+b(n)*sin(n*x);
    f=f+a(n)*cos(pi/2*n*x)+b(n)*sin(omega(n)*x);
% f=f+a(n)*cos(pi*n*x)+b(n)*sin(pi*n*x);
end
close all;
% figure("Position",[100,100,297,210])
% plot(x,f); hold on;
% % plot(x,exp(x));
% plot(x,x/4);
plot(x,f);
% ylim([min(f) max(f)])
%  xlim([min(0) max(x)]); 
grid on;
%%
% exportgraphics(gcf,"plot"+num2str(nm)+".png","Resolution",300);
exportgraphics(gcf,"plot5-"+num2str(nm)+".png","Resolution",300);