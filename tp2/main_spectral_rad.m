%% tum
A=[0.7, -0.4; -0.2, 0.5];
D=diag(diag(A));
% eig(eye(2)-inv(D)*A)
eig(eye(2)-A)
%  eig(eye(2)-inv(+diag(1,-1)*diag(A,-1)+D)*A)
%  om=1.08; eig(eye(2)-inv(diag(1,-1)*diag(A,-1)+D/om)*A)

M=eye(2)-inv(D)*A %jacobi
eig(M)'
om=2/(2-max(eig(M))-min(eig(M)))
eig(eye(2)-inv(D/om)*A)'

%% su
A=[1 2 3; 2 5 10; 3 10 26];
D=diag(diag(A));

om=2/(2-max(eig(A))-min(eig(A)));
eig(eye(3)-inv(D)*A)
% eig(eye(3)-A)
 
 eig(eye(3)-inv(+diag(diag(A,-1),-1)+D)*A)
 om=0.5; eig(eye(3)-inv(+diag(diag(A,-1),-1)+D/om)*A)
%% su 2
A=[4, -1,0; -1,4,-1;0,-2,4];
D=diag(diag(A));

eig(eye(3)-inv(D)*A)'
om=2/(2-max(eig(A))-min(eig(A)))
eig(eye(3)-inv(D/om)*A)'

% eig(eye(3)-A)
  evs=eig(eye(3)-inv(+diag(diag(A,-1),-1)+D)*A)
  rho=max(max(evs));
  om=2/(1+sqrt(1-rho^2))
 om=1.05;
 eig(eye(3)-inv(+diag(diag(A,-1),-1)+D/om)*A)

% M=eye(3)-inv(D)*A %jacobi
% eig(M)'
% om=2/(2-max(eig(M))-min(eig(M)))
% eig(eye(3)-inv(D/om)*A)'

 %%
 hold on; axis equal; grid on;
 for om=0.9:0.01:1.1
     v=eig(eye(3)-inv(+diag(diag(A,-1),-1)+D/om)*A);
     plot(real(v),imag(v)+om,'x')
 end
