function [x,h,xw,xe] = mesh_properties(L,Nx)
% MESH_PROPERTIES Computes geometrical properties required for the finite volume method. 

x = linspace(0,L,Nx)'; % Node locations
h = L/(Nx-1);
dx = diff(x); % Node spacing
xw = zeros(Nx,1); % West control volume boundary locations
xe = zeros(Nx,1); % East control volume boundary locations

xw(1) = 0; % Left boundary node
xe(1) = dx(1)/2; 
for i = 2:Nx-1 % Internal nodes
    xw(i) = x(i) - dx(i-1)/2;
    xe(i) = x(i) + dx(i)/2;
end
xw(Nx) = x(Nx) - dx(Nx-1)/2; % Right boundary node
xe(Nx) = x(Nx);