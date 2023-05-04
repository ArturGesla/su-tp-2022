clc; clear;
Q=[];
Qm=[];
fluxVarray=[];
%%
% nx=23;
nx=47;
% nx=95;
% nx=191;
ny=(nx+1)/1.2-1;
nn=nx*ny;
h=1/(nx+1);
A=zeros(nn);
b=zeros(nn,1);
omega=1.8;
gs = 0.; gw = 1.; ge = 0.2; gn = 0.6; f = 0;

for i=1:nn
    A(i,i)=4;
    if (i<nn && mod(i,nx)~=0) A(i,i+1)=-1; end
    if (i>1 && mod(i-1,nx)~=0) A(i,i-1)=-1; end
    if(i>nx)A(i,i-nx)=-1; end
    if(i<(ny-1)*nx+1)A(i,i+nx)=-1; end
    if(mod(i-1,nx)==0)b(i)=b(i)+gw; end
    if(i<nx+1)b(i)=b(i)+gs; end
    if(mod(i,nx)==0)b(i)=b(i)+ge; end
    if(i>(ny-1)*nx)b(i)=b(i)+gn; end
end
% A=sparse(A);
q=ones(nn,1);
r=1;

% % for outeri=1:175
% outeri=1;
% while (norm(r/nn)>1e-5)
%     outeri=outeri+1;
% r=b-A*q;
% dq=zeros(nn,1);
% for i=1:nn
%     dq(i)=omega*r(i);
%     for j=1:i-1
%  dq(i)=dq(i)-A(i,j)*omega*dq(j);
%     end
%     dq(i)=dq(i)/A(i,i);
% end
% q=q+dq;
% fprintf("%d\t%4.2e\n",outeri,norm(r/nn))
% end

%
q1=reshape(q,[nx,ny]);
[xc,zc]=meshgrid([1:nx]./nx,[1:ny]./ny);
x=sparse(A)\b;
% x=q;
q2=reshape(x,[nx,ny])';
iy=(ny+1)/2
dtdy=[0,(q2(iy+1,:)-q2(iy-1,:))/2/h, 0]
fluxV=trapz(dtdy)*h
fluxVarray=[fluxVarray,fluxV]
% Q=[Q,norm(x)/nx]
% Qm=[Qm,x((nx-1)/2*nx+(nx+1)/2)]
%%
% (Q(end-1)-Q(end-2))/(Q(end)-Q(end-1))
% (Qm(end-1)-Qm(end-2))/(Qm(end)-Qm(end-1))
(fluxVarray(end-1)-fluxVarray(end-2))/(fluxVarray(end)-fluxVarray(end-1))