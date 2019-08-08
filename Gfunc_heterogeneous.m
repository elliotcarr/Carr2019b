function G = Gfunc_heterogeneous(t,u,D,h,xe,xw,Nx,g0,gL)
% GFUNC_HETEROGENEOUS Right-hand side function passed to ode15s to solve
% the heterogeneous model (6)-(8) using a finite volume discretisation.

G = zeros(Nx, 1);

G(1) = u(1) - g0(t); % Left boundary node
for i = 2:Nx-1 % Interior nodes
    G(i) = D(xe(i))*(u(i+1) - u(i))/h^2 - D(xw(i))*(u(i) - u(i-1))/h^2;
end
G(Nx) = u(Nx) - gL(t); % Right boundary node
