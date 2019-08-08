close all; clc; clear all;

Nx = 1001; % Number of nodes to solve BVPs and homogenized/heterogeneous models

% paper_formatting = true;
paper_formatting = false;

font_size = 34;
if paper_formatting
    path_name = '../../Paper/Figures/';
    addpath('../export_fig-master')
end

%% Specify problem
L = 1; % Length of medium
D = @(x) 0.8 - 0.6*x + 0.2*sin(20*x); % Spatially-varying (heterogeneous) diffusivity
f = @(x) zeros(size(x)); % Initial condition
g0 = @(t) 1*ones(size(t)); % Boundary value at x = 0
gL = @(t) 0*ones(size(t)); % Boundary value at x = L
tspan = [1e-2,1e-1,1e0]; % Plot solution at these times
Nt = 100; terror = linspace(0,1,Nt+1); terror(1) = []; % Calculate error based on these times

% Test cases from Carr (2019)
% run('./Test Cases/CaseA');
run('./Test Cases/CaseB');
% run('./Test Cases/CaseC');
% run('./Test Cases/CaseD');
% run('./Test Cases/CaseE');
% run('./Test Cases/CaseF');

%% Finite volume geometrical properties
[x,h,xw,xe] = mesh_properties(L,Nx);

%% Plot diffusivity D(x)
if paper_formatting
    figure;
    set(gcf,'Color','w')    
    if strcmp(CaseNum,'E')
        for i = 2:length(Dvec)-1
            plot([(i-2)/16,(i-1)/16],[Dvec(i),Dvec(i)],'LineWidth',2,'Color',[0,0.447,0.741])
            hold on            
            plot([(i-1)/16,(i-1)/16],[Dvec(i),Dvec(i+1)],'--','LineWidth',2,'Color',[0,0.447,0.741])
        end
    else
        plot(x,D(x),'LineWidth',2);
    end
    axis([0 1 0 1])
    set(gca,'Fontsize',font_size-2,'TickLabelInterpreter','latex','Xtick',[0,1],'YTick',[0,1])
    xlabel('$x$','Interpreter','LaTeX','Fontsize',font_size)
    ylabel('$D(x)$','Interpreter','LaTeX','Fontsize',font_size)
    if ~(strcmp(CaseNum,'Intro1') || strcmp(CaseNum,'Intro2'))
        title([FigLabels{1},' Case ',CaseNum],'Interpreter','LaTeX','Fontsize',font_size)
    end
    vec_pos = get(get(gca, 'XLabel'), 'Position');
    set(get(gca, 'XLabel'), 'Position', vec_pos + [0 0.15 0]);
    vec_pos = get(get(gca, 'YLabel'), 'Position');
    set(get(gca, 'YLabel'), 'Position', vec_pos + [0.05 0 0]);
    drawnow
    feval('export_fig',[path_name,'diffusivity_',CaseNum],'-pdf')
end

%% Solve boundary value problems for s(x) and w(x) over heterogeneous medium
r = zeros(size(x));
s0 = g0(Inf);
sL = gL(Inf);
s = bvp_heterogeneous(D,h,xw,xe,Nx,r,s0,sL);

fx = f(x);
r = fx - s;
w0 = integral(@(t) g0(Inf)-g0(t),0,Inf);
wL = integral(@(t) gL(Inf)-gL(t),0,Inf);
w = bvp_heterogeneous(D,h,xw,xe,Nx,r,w0,wL);

%% Define nonlinear system
S0 = s0; SL = sL; W0 = w0; WL = wL; % Homogenized boundary values match those imposed on heterogeneous model
wint = trap_rule(h,w);
sint = trap_rule(h,s);
F = @(coeffs) Ffunc(coeffs(1),coeffs(2),Nx,h,fx,sint,wint,S0,SL,W0,WL);

%% Solve nonlinear system for effective coefficients
r = zeros(size(x));
Deff_tilde = L / integral(@(x) 1./D(x),0,L);
options = optimoptions('fsolve','TolX',1e-20,'TolFun',1e-24,'Display','none');
[eff,fval] = fsolve(F,[Deff_tilde,0.1*Deff_tilde],options);
Ffinal = max(abs(F(eff)));
Deff = eff(1);
veff = eff(2);

%% Solve homogenized model with effective parameters
M = eye(Nx); M(1,1) = 0; M(Nx,Nx) = 0; options = odeset('Mass',M);
[~,U] = ode15s(@(t,U) Gfunc_homogenized(t,U,Deff,veff,h,Nx,g0,gL),[0,tspan],fx,options);
[~,Utilde] = ode15s(@(t,U) Gfunc_homogenized(t,U,Deff_tilde,0.0,h,Nx,g0,gL),[0,tspan],fx,options);

%% Solve heterogeneous model
[~,u] = ode15s(@(t, u) Gfunc_heterogeneous(t,u,D,h,xe,xw,Nx,g0,gL),[0,tspan],fx,options);

%% Plot solutions of homogenized and heterogeneous models
figure;
set(gcf,'Color','w')
for i = 2:length(tspan)+1 % Don't plot initial condition
    p1 = plot(x,u(i,:),'-','Color','k','LineWidth',1);
    hold on
    p2 = plot(x,Utilde(i,:),'--','Color','k','LineWidth',1);
    p3 = plot(x,U(i,:),'-','Color',[204,37,41]/255,'LineWidth',3);
end
axis([0 1 0 1])
set(gca,'Fontsize',font_size-2,'TickLabelInterpreter','latex','Xtick',[0,1],'YTick',[0,1])
xlabel('$x$','Interpreter','LaTeX','Fontsize',font_size)
ylabel('Solution','Interpreter','LaTeX','Fontsize',font_size)
if paper_formatting
    title([FigLabels{2},' Case ',CaseNum],'Interpreter','LaTeX','Fontsize',font_size)
else
    leg = legend([p1(1),p2(1),p3(1)],'$u(x,t)$','$\widetilde{U}(x,t)$','$U(x,t)$','Location',...
        'NorthEast','Orientation','vertical');
    set(leg,'Fontsize',font_size-6,'Interpreter','LaTeX')
end
vec_pos = get(get(gca, 'XLabel'), 'Position');
set(get(gca, 'XLabel'), 'Position', vec_pos + [0 0.15 0]);
vec_pos = get(get(gca, 'YLabel'), 'Position');
set(get(gca, 'YLabel'), 'Position', vec_pos + [0.05 0 0]);
drawnow
if paper_formatting
    feval('export_fig',[path_name,'solution_homogenized_',CaseNum],'-pdf')
end

%% Compute errors
[~,u] = ode15s(@(t, u) Gfunc_heterogeneous(t,u,D,h,xe,xw,Nx,g0,gL),[0,terror],fx,options);
[~,U] = ode15s(@(t,U) Gfunc_homogenized(t,U,Deff,veff,h,Nx,g0,gL),[0,terror],fx,options);
[t,Utilde] = ode15s(@(t,U) Gfunc_homogenized(t,U,Deff_tilde,0.0,h,Nx,g0,gL),[0,terror],fx,options);

epsilon = 0; % homogenized advection-diffusion model
epsilon_tilde = 0; % homogenized diffusion-only model
for j = 2:Nt+1 % Don't include initial solution
    epsilon = epsilon + sum(abs(u(j,:)-U(j,:)));
    epsilon_tilde = epsilon_tilde + sum(abs(u(j,:)-Utilde(j,:)));    
end
epsilon = epsilon/(Nx*Nt);
epsilon_tilde = epsilon_tilde/(Nx*Nt);

if paper_formatting
    fprintf('\\num{%1.3f} & \\num{%1.3f} & \\num{%1.2e} & \\num{%1.3f} & \\num{%1.2e}\n',...
        Deff,veff,epsilon,Deff_tilde,epsilon_tilde)
end