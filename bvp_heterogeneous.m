function m = bvp_heterogeneous(D,h,xw,xe,Nx,r,m0,mL)
% BVP_HETEROGENEOUS Solves one-dimensional steady state heterogeneous
% diffusion problems (13)-(14) and (20)-(22) using a finite volume
% discretisation.

A = zeros(Nx,Nx);
b = h*r;

A(1,1) = 1; % Left boundary node
b(1) = m0;
for i = 2:Nx-1  % Interior nodes
    A(i,i-1) = D(xw(i))/h;
    A(i,i) = -D(xw(i))/h - D(xe(i))/h;
    A(i,i+1) = D(xe(i))/h;
end
A(Nx,Nx) = 1; % Right boundary node
b(Nx) = mL;

m = A\b;