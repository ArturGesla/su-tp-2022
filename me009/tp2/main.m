nx=3;
ny=2;
nn=nx*ny;
A=zeros(nn);
b=zeros(nn,1);
omega=1.9;
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
A
q=ones(nn,1)
r=b-A*q

dq=zeros(nn,1);
for i=1:nn
    dq(i)=omega*r(i);
    for j=1:i-1
 dq(i)=dq(i)-A(i,j)*omega*dq(j);
    end
    dq(i)=dq(i)/A(i,i);
end