

a=@(n)(-1)^n*(exp(pi)-exp(-pi))/pi/(n^2+1);
b=@(n)-n*a(n);

x=[-2*pi:0.01:2*pi];
f=a(0)/2;
nm=50;
for n=1:nm
    f=f+a(n)*cos(n*x)+b(n)*sin(n*x);
end
figure("Position",[100,100,297,210])
plot(x,f); hold on;
plot(x,exp(x));
ylim([min(f) max(f)])
xlim([min(x) max(x)]); grid on;
exportgraphics(gcf,"plot"+num2str(nm)+".png","Resolution",300);