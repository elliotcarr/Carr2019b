% Problem parameters for Case D

L = 1; % Length of medium
D = @(x) 0.1*exp(2.2*x.^2) + 0.05*sin(80*x); % Spatially-varying (heterogeneous) diffusivity
f = @(x) zeros(size(x)); % Initial condition
g0 = @(t) 1*ones(size(t)); % Boundary value at x = 0
gL = @(t) 0*ones(size(t)); % Boundary value at x = L
tspan = [1e-2,1e-1,1e0]; % Plot solution at these times
Nt = 100; terror = linspace(0,1,Nt+1); terror(1) = []; % Calculate error based on these times
CaseNum = 'D'; FigLabels = {'(d)','(m)'};