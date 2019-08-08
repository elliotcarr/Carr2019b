function F = Ffunc(Deff,veff,Nx,h,fx,sint,wint,S0,SL,W0,WL)
% FFUNC Vector-valued function defining the nonlinear system (26)-(28).

F = zeros(2,1);

r = zeros(Nx,1);
S = bvp_homogenized(Deff,veff,h,Nx,r,S0,SL);
F(1) = trap_rule(h,S) - sint;

r = fx - S;
W = bvp_homogenized(Deff,veff,h,Nx,r,W0,WL);
F(2) = trap_rule(h,W) - wint;