function I = trap_rule(h,f)
% TRAP_RULE Trapezoidal rule approximation. 

I = h/2*(sum(f(1:end-1)+f(2:end)));