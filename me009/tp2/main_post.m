t=[]; h=[]; nm=[];
%%
a=importdata("isoP_sor.dat")'; n=length(a);%nx=30; ny=25;
% t=[t,a((n+1)/2,(n+1)/2)]
% h=[h,1.5/(n+1)]
% nm=[nm,norm(a,"fro")*1.5/(n+1)]
contour(a)
%%
(nm(end-1)-nm(end-2))/((nm(end)-nm(end-1)))
%%
q0=-4.58566070e-01;
q1=-4.56434309e-01
q2=-4.55836892e-01
q3=-4.55676645e-01
% q4=-4.55582261e-01
(q1-q0)/(q2-q1)
(q2-q1)/(q3-q2)
% (q3-q2)/(q4-q3)
% plot([q1,q2,q3,q4])

%% analytical sol
close all;
figure("Position",[100,100,1200,600])
l=1.25;
L=1.5;

b0=1/2*l/L;
% bn=@(k)2/l/sinh(k*pi)*l^2/k/k/pi/pi*(-1+(-1)^k);
bn=@(k)2/l*l^2/k/k/pi/pi*(-1+(-1)^k);

% nx=5; ny=4;
nx=23; ny=19;
% nx=47; ny=39;
% nx=95; ny=79;
[x,y]=meshgrid(linspace(0,L,nx+2),linspace(0,l,ny+2));
% x=1.25;
% y=x;
T=b0*x;
nmodes=100;
for i=1:nmodes
    T=T+bn(i).*cos(y*i*pi/l).*sinh(i*pi/l*x)./sinh(i*pi*L/l);
end
subplot(1,2,1);
contour(x,y,T,20); colorbar();
contour(x(:,2:end-1),y(:,2:end-1),T(:,2:end-1),20); colorbar();

a=importdata("isoP_sor.dat")'; n=length(a);%nx=30; ny=25;
% t=[t,a((n+1)/2,(n+1)/2)]
% h=[h,1.5/(n+1)]
% nm=[nm,norm(a,"fro")*1.5/(n+1)]
subplot(1,2,2)
contour(x(:,2:end-1),y(:,2:end-1),a,20); colorbar();
close all;
mesh(T(:,2:end-1)-a)
% mesh(t1-T)
% mesh(T)

%% etu
a=importdata("isoP_sor.dat");
x=reshape(a(:,1),[97,81])
y=reshape(a(:,2),[97,81])
t2=reshape(a(:,3),[97,81])
mesh(x,y,t2)

T=b0*x;
nmodes=50;
for i=1:nmodes
    T=T+bn(i).*cos(y*i*pi/l).*sinh(i*pi/L*x);
end
mesh(x,y,T-t2)