% Problem parameters for Case E

L = 1; % Length of medium
load CaseE.mat
D = griddedInterpolant(xvec,Dvec,'nearest'); % Spatially-varying (heterogeneous) diffusivity
f = @(x) exp(-30*(x-0.5).^2); % Initial condition
g0 = @(t) 0*ones(size(t)); % Boundary value at x = 0
gL = @(t) 1e-6*ones(size(t)); % Boundary value at x = L
tspan = [1e-2,1e-1,1e0]; % Plot solution at these times
Nt = 100; terror = linspace(0,1,Nt+1); terror(1) = []; % Calculate error based on these times
CaseNum = 'E'; FigLabels = {'(e)','(k)'};