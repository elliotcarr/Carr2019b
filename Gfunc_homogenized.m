function G = Gfunc_homogenized(t,U,Deff,veff,h,Nx,g0,gL)
% GFUNC_HOMOGENIZED Right-hand side function passed to ode15s to solve
% the homogenized model (9)-(11) using a finite volume discretisation.

G = zeros(Nx,1);

G(1) = U(1) - g0(t); % Left boundary node
for i = 2:Nx-1 % Interior nodes
    G(i) = Deff*(U(i+1) - U(i))/h^2 - Deff*(U(i) - U(i-1))/h^2 + veff*(U(i-1) - U(i+1))/(2*h);
end
G(Nx) = U(Nx) - gL(t); % Right boundary node
