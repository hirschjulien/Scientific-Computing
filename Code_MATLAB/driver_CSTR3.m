% driver_CSTR3 Script
%
%
clc
clear
% --------------------
% Parameters
% --------------------
% p.Nz = 100;                     % [-]           Resolution of C in the spatial dimension
p.rho  = 1.0;                   % [kg/L]        Density
p.c_p  = 4.186;                 % [kj/(kg*K)]   Specific heat capacity
p.k0   = exp(24.6);             % [L/(mol*s)]   Arrhenius constant
p.EaR  = 8500.0;                % [K]           Activation energy
p.DHr  = -560.0;                % [kJ/mol]      Enthalpy of reaction
p.L    = 10;                    % [m]           Length of reactor
p.A    = 0.1;                   % [m^2]         Cross-sectional area of reactor
p.V    = 0.105; %p.A*p.L*1000;          % [L]           Reactor volume (A*L=1m^3)(1m^3=1000L)
p.DA   = [0.1; 0.1; 0.1];       % [m^2/s]       Diffusion coefficients
p.beta = -p.DHr/(p.rho*p.c_p);
p.nu = [-1; -2; p.beta];        % 3 state only!
% p.dz = p.L/p.Nz;                % PFR case only!
% Boundary conditions
cAin = 1.6/2;                   % [mol/L]       Concentration of A at Inlet
cBin = 2.4/2;                   % [mol/L]       Concentration of B at Inlet
Tin  = 273.65;                  % [K]           Temperature at Inlet (=0.5C)
u=[cAin;cBin;Tin];              % 3 state case
p.u=u;
% --------------------
 p.F=650;                       % [mL/min]     Volumetric Flowrate

% Create the parameter struct
%v = F/A not p.v;
% u = Tin;                      % Boundary condition
InitCond = [0;0;Tin];           % Initial conditions
%InitCond = [0.5;1;250]; 
t0 = 0;
tf = 1;                         %

% aa = linspace(0, 1000, 2500); % X-axis
% mu = 500;                     % Mean
% sigma = 100;                  % Standard deviation
% p.F = (1/(sigma * sqrt(2*pi))) * exp(-0.5 * ((aa - mu)/sigma).^2);
set(groot, 'defaultLineLineWidth', 1.8);
%% ODE45
options = odeset('RelTol',1.0e-7,'AbsTol',1.0e-7);
[T,X]= ode45(@CSTR3, [t0, tf], InitCond, options, u, p);


% MAKE THE plotter_CSTR.m
%% PLOTS FOR ODE45
figs1 = tiledlayout(3,1);
figs1.Padding = 'compact';
figs1.TileSpacing = 'compact';
ax1 = nexttile;
plot(T,X(:,1));
%pbaspect(ax1, [1 1 1])
title('Solved with ode45');
ylim([0,1]);
ylabel('$$C_A(t)$$','interpreter','latex')
ax2 = nexttile;
plot(T,X(:,2));
ylim([0,1]);
ylabel('$$C_B(t)$$','interpreter','latex');
ax3 = nexttile;
plot(T,X(:,3));
hold on
plot(t1,x1);
legend("CSTR3", "CSTR1");
ylabel('$$T(t)$$','interpreter','latex');
xlabel('$$t [min]$$', 'interpreter', 'latex') 
exportgraphics(figs1,'CSTR_ode45_constflow.pdf','ContentType','vector')
disp("Fin!")

%% ode15s
clc % Clears all text from command window
options = odeset('Jacobian',@CSTR3Jac,'RelTol',1.0e-14,'AbsTol',1.0e-14);
[T,X]=ode15s(@CSTR3, [t0, tf], InitCond, options, p);

%% Plots for ode15s
figs3 = tiledlayout(3,1);
figs3.Padding = 'compact';
figs3.TileSpacing = 'compact';
ax4 = nexttile;
plot(T,X(:,1));
%pbaspect(ax1, [1 1 1])
title('Solved with ode15s');
xlim([0, 0.001]);
ylim([0,1]);
ylabel('$$C_A(t)$$','interpreter','latex')
ax5 = nexttile;
plot(T,X(:,2));
ylim([0,1]);
ylabel('$$C_B(t)$$','interpreter','latex');
ax6 = nexttile;
plot(T,X(:,3));
hold on
plot(Tcstr1,Xcstr1);
legend("CSTR3", "CSTR1");
xlim([0, 0.001]);
ylabel('$$T(t)$$','interpreter','latex');
xlabel('$$t [min]$$', 'interpreter', 'latex') 
exportgraphics(figs3,'CSTR_ode15s_constflow.pdf','ContentType','vector')
disp("Fin!")

%% ImplicitEulerFixed
clc
N = 10^4;            % Number of steps, with 10^3 steps it does not converge.   
%options = odeset('RelTol',1.0e-6,'AbsTol',1.0e-6); % The options are not needed as input in the custom ExplicitEuler function
[T,X]=ImplicitEulerFixedStepSize(@CSTR3FunJac, [t0, tf], N, 1.e-6, InitCond, p);
%%
figs3 = tiledlayout(3,1);
figs3.Padding = 'compact';
figs3.TileSpacing = 'compact';
ax31 = nexttile;
plot(T,X(:,1));
%pbaspect(ax1, [1 1 1])
title('Implicit Euler fixed step');
ylim([0,1]);
ylabel('$$C_A(t)$$','interpreter','latex')
ax32 = nexttile;
plot(T,X(:,2));
ylim([0,1]);
ylabel('$$C_B(t)$$','interpreter','latex');
ax33 = nexttile;
plot(T,X(:,3));
hold on
plot(Tcstr1,Xcstr1);
legend("CSTR3", "CSTR1");
ylabel('$$T(t)$$','interpreter','latex');
xlabel('$$t [min]$$', 'interpreter', 'latex') 
% exportgraphics(figs3,'CSTR_impleul_constflow.pdf','ContentType','vector')
% disp("Fin!")

%% ImplicitEulerAdaptive
clc
abstol = 1e-6; % Usual / default value 
reltol = 1e-6;
% INITIAL STEP SIZE CALCULATION
h0 = InitialStepSize(@CSTR3,t0,InitCond,p);

[T,X]=ImplicitEulerAdaptiveStep(@CSTR3FunJac, [t0, tf], InitCond, h0, abstol, reltol, p);

figs3 = tiledlayout(3,1);
figs3.Padding = 'compact';
figs3.TileSpacing = 'compact';
ax31 = nexttile;
plot(T,X(:,1));
%pbaspect(ax1, [1 1 1])
title('Implicit Euler adaptive step');
ylim([0,1]);
ylabel('$$C_A(t)$$','interpreter','latex')
ax32 = nexttile;
plot(T,X(:,2));
ylim([0,1]);
ylabel('$$C_B(t)$$','interpreter','latex');
ax33 = nexttile;
plot(T,X(:,3));
hold on
plot(Tcstr1,Xcstr1);
legend("CSTR3", "CSTR1");
ylabel('$$T(t)$$','interpreter','latex');
xlabel('$$t [min]$$', 'interpreter', 'latex') 
% exportgraphics(figs3,'CSTR_impleuladaptve_constflow.pdf','ContentType','vector')
% disp("Fin!")

%%  RK4 Fixed, ClasicalRungeKuttaFixedStep.m
N=1e6;

[T,X]=ClassicalRungeKuttaFixedStep(@CSTR3, [t0, tf], InitCond, N, p);

figs3 = tiledlayout(3,1);
figs3.Padding = 'compact';
figs3.TileSpacing = 'compact';
ax31 = nexttile;
plot(T,X(:,1));
%pbaspect(ax1, [1 1 1])
title('RK4 fixed step');
ylim([0,1]);
ylabel('$$C_A(t)$$','interpreter','latex')
ax32 = nexttile;
plot(T,X(:,2));
ylim([0,1]);
ylabel('$$C_B(t)$$','interpreter','latex');
ax33 = nexttile;
plot(T,X(:,3));
hold on
plot(Tcstr1,Xcstr1);
legend("CSTR3", "CSTR1");
ylabel('$$T(t)$$','interpreter','latex');
xlabel('$$t [min]$$', 'interpreter', 'latex') 
% exportgraphics(figs3,'CSTR_RK4fixed_constflow.pdf','ContentType','vector')
% disp("Fin!")

%% ADAPTIVE, CLassicalRungeKuttaAdaptiveStep.m +
% VAriables needed for adaptive step method:
abstol = 1e-4; % odeset default value is 1e-6
reltol = 1e-4; % odeset default is 1e-3
% INITIAL STEP SIZE CALCULATION
h0 = InitialStepSize(@CSTR3,t0,InitCond,p);

[T,X]=ClassicalRungeKuttaAdaptiveStep(@CSTR3, [t0, tf], InitCond, h0, abstol, reltol, p);

figs3 = tiledlayout(3,1);
figs3.Padding = 'compact';
figs3.TileSpacing = 'compact';

ax31 = nexttile;
plot(T,X(:,1));
%pbaspect(ax1, [1 1 1])
% xlim([0 0.001])
% ylim([0,1]);
title('RK4 adaptive step');
ylabel('$$C_A(t)$$','interpreter','latex')
ax32 = nexttile;
plot(T,X(:,2));
% xlim([0 0.001])
% ylim([0,1]);
ylabel('$$C_B(t)$$','interpreter','latex');
ax33 = nexttile;
plot(T,X(:,3));
% xlim([0 0.001])
% hold on
% plot(Tcstr1,Xcstr1);
legend("CSTR3", "CSTR1");
ylabel('$$T(t)$$','interpreter','latex');
xlabel('$$t [min]$$', 'interpreter', 'latex') 
% exportgraphics(figs3,'CSTR_RK4adaptve_constflow.pdf','ContentType','vector')
% disp("Fin!")


%% Many flow rates ADAPTIVE, CLassicalRungeKuttaAdaptiveStep.m
% Variables needed for adaptive step method:
abstol = 1e-6; % Usual / default value is 1e-6
reltol = 1e-6; % Default is 1e-3 in odeset.
% INITIAL STEP SIZE CALCULATION
h0 = InitialStepSize(@CSTR1,t0,InitCond,p);
% h0=0.0001;

Ff(1) = 100;    % Creating an array with the flow rates
Ff(2) = 250;
Ff(3) = 400;
Ff(4) = 550;
Ff(5) = 700;
Ff(6) = 850;
Ff(7) = 1000;
Tcstr3=cell(1,7);
Xcstr3=cell(1,7);

for k=1:7
    p.F=Ff(k);
    [Tcstr3{k},Xcstr3{k}]=ClassicalRungeKuttaAdaptiveStep(@CSTR3, [t0, tf], InitCond, h0, abstol, reltol, p);
end

%% Plots for rk4 adaptive many constant flows comparison
figure('Position', [100 100 1000 500]);  % [left, bottom, width, height]
figs3 = tiledlayout(3,1);
figs3.Padding = 'compact';
figs3.TileSpacing = 'compact';

ax41 = nexttile;
for k=1:7
    plot(Tcstr3{k}, Xcstr3{k}(:,1) )
    hold on
end
grid on
pbaspect(ax41, [5 1 1])
xlim([0 0.01])
ylim([0,1]);
title('RK4 adaptive step');
ylabel('$$C_A(t)$$','interpreter','latex')

ax42 = nexttile;
for k=1:7
    plot(Tcstr3{k}, Xcstr3{k}(:,2) )
    hold on
end
grid on
pbaspect(ax42, [5 1 1])
xlim([0 0.01])
ylim([0,1]);
ylabel('$$C_B(t)$$','interpreter','latex');

ax43 = nexttile;
for k=1:7
    plot(Tcstr3{k}, Xcstr3{k}(:,3), "DisplayName", "F="+string(Ff(k)) )
    hold on
end
grid on
pbaspect(ax43, [5 1 1])
xlim([0 0.01])

% yticks = linspace(min(Xcstr3{1}(:,3)), max(Xcstr3{1}(:,3)), 4);
% ax43.YTick = yticks;
% ax43.YTickLabel = sprintfc('%.1e', yticks - 273.65);

lgd=legend("location", "eastoutside");
lgd.Title.String = '[mL/min]';d
ylabel('$$T(t)$$','interpreter','latex');
xlabel('$$t [min]$$', 'interpreter', 'latex') 
exportgraphics(figs3,'cstr3_RK4_adapt_flows.pdf','ContentType','vector')
disp("Fin!")



%%
for k=1:7
    plot(Tcstr3{k}, Xcstr3{k}, "DisplayName", "F="+string(Ff(k)) )
    hold on
end
hold off
pbaspect([5 1 1]);
lgd=legend("location", "eastoutside");
lgd.Title.String = '[mL/min]';

xlim([0,0.01])

ax=gca;
yticks = linspace(min(Xcstr3{1}), max(Xcstr3{1}), 4);
ax.YTick = yticks;
ax.YTickLabel = sprintfc('%.1e', yticks - 273.65);

title("CSTR3, RK4 Adaptive")
xlabel('$$t[min]$$','interpreter','latex')
ylabel('$$T[C]$$','interpreter','latex')
grid on
% exportgraphics(fig3B,'cstr1_RK4_adapt_flows.pdf','ContentType','vector') % X in "XCV" is the factor multiplied to the Correct Volume used
% disp('graph exported')