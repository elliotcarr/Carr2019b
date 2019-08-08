function M = bvp_homogenized(Deff,veff,h,Nx,R,M0,ML)
% BVP_HETEROGENEOUS Solves one-dimensional steady state homogenized
% diffusion problems (15)-(16) and (23)-(25) using a finite volume
% discretisation.

A = zeros(Nx,Nx);
b = h.*R;

A(1,1) = 1; % Left boundary node
b(1) = M0;
for i = 2:Nx-1 % Interior nodes
    A(i,i-1) = Deff/h + veff/2;
    A(i,i) = -Deff/h - Deff/h;
    A(i,i+1) = Deff/h - veff/2;
end
A(Nx,Nx) = 1; % Right boundary node
b(Nx) = ML; 

M = A\b;