% Problem parameters for Case F

L = 1; % Length of medium
load CaseF.mat
D = griddedInterpolant(xvec,Dvec,'linear'); % Spatially-varying (heterogeneous) diffusivity
f = @(x) (x<0.5).*(2*x) + (x>=0.5).*(2-2*x); % Initial condition
g0 = @(t) 0*ones(size(t)); % Boundary value at x = 0
gL = @(t) 0.5*ones(size(t)); % Boundary value at x = L
tspan = [1e-2,1e-1,1e0]; % Plot solution at these times
Nt = 100; terror = linspace(0,1,Nt+1); terror(1) = []; % Calculate error based on these times
CaseNum = 'F'; FigLabels = {'(f)','(l)'};