% driver_CSTR1 Script
%
%
clc
clear
% --------------------
% Parameters
% --------------------
p.Nz = 100;                     % [-]           Resolution of C in the spatial dimension
p.rho  = 1.0;                   % [kg/L]        Density
p.c_p  = 4.186;                 % [kj/(kg*K)]   Specific heat capacity
p.k0   = exp(24.6);             % [L/(mol*s)]   Arrhenius constant
p.EaR  = 8500.0;                % [K]           Activation energy
p.DHr  = -560.0;                % [kJ/mol]      Enthalpy of reaction
p.V    = 0.105;          % [L]           Reactor volume (A*L=1m^3)(1m^3=1000L)
p.DA   = [0.1; 0.1; 0.1];       % [m^2/s]       Diffusion coefficients
p.beta = -p.DHr/(p.rho*p.c_p);
% p.nu = [-1; -2; p.beta];        % 3 state only!
% p.dz = p.L/p.Nz;                % PFR case only!
% Boundary conditions
cAin = 1.6/2;                   % [mol/L]       Concentration of A at Inlet
cBin = 2.4/2;                   % [mol/L]       Concentration of B at Inlet
Tin  = 273.65;                  % [K]           Temperature at Inlet (=0.5C)
u=[cAin;cBin;Tin];              % 3 state case
p.u=u;
% --------------------
p.F=650;                        % [mL/min]     Volumetric Flowrate

% u = Tin;      % Boundary condition
T0 = Tin;       % Initial conditions
% Time interval
t0 = 0;
tf = 1;      % 35*60 :Data point every 10 seconds

%%
options = odeset('RelTol',1.0e-10,'AbsTol',1.0e-12);
[T1,X1]= ode45(@CSTR1, [t0, tf], T0, options, p);

fig1= figure(1);
plot(T1,X1)
xlabel('$$t$$','interpreter','latex')
ylabel("Temperature")
grid on

%% Implicit solver - ode15s
clc
options = odeset('Jacobian',@CSTR1Jac,'RelTol',1.0e-10,'AbsTol',1.0e-12);
[Tcstr1,Xcstr1]=ode15s(@CSTR1, [t0, tf], T0, options, p);

fig2= figure(2);
plot (Tcstr1,Xcstr1)
grid on

%% ImplicitEulerFixedStep
clc
N = 10^4;            % Number of steps, with 10^3 steps it does not converge.   
%options = odeset('RelTol',1.0e-6,'AbsTol',1.0e-6); % The options are not needed as input in the custom ExplicitEuler function
[Tcstr1,Xcstr1]=ImplicitEulerFixedStepSize(@CSTR1FunJac, [t0, tf], N, 1.e-6, T0, p);

fig3= figure(1);
plot(Tcstr1,Xcstr1)
xlabel('$$t[min]$$','interpreter','latex')
ylabel("Temperature")
grid on

%% ADAPTIVE ImplicitEuelrAdaptiveStep.m + 
clc
abstol = 1e-10; % Usual / default value 
reltol = 1e-12;
% INITIAL STEP SIZE CALCULATION
h0 = InitialStepSize(@CSTR1,t0,T0,p);

[Tcstr1,Xcstr1]=ImplicitEulerAdaptiveStep(@CSTR1FunJac, [t0, tf], T0, h0, abstol, reltol, p);

%%
figs4 = tiledlayout(2,1);
figs4.Padding = 'compact';
figs4.TileSpacing = 'compact';
ax31 = nexttile;
plot(Tcstr1,Xcstr1)
pbaspect([5, 1, 1])
xlim([0,2100])
grid on
ylabel('$$T[K]$$','interpreter','latex')
title("Impl.Euler adaptive, atol="+abstol+", rtol="+reltol)
ax31 = nexttile;
plot(Tcstr1,Xcstr1)
pbaspect([5, 1, 1])
xlim([0,.01])
grid on
xlabel('$$t[s]$$','interpreter','latex')
ylabel('$$T[K]$$','interpreter','latex')
% exportgraphics(figs3,'cstr1_CV_impleuladapt.pdf','ContentType','vector')
% % X in "XCV" is the factor multiplied to the Correct Volume used
% disp('graph exported')
%% RK4 Fixed, ClasicalRungeKuttaFixedStep.m

N=1e7;

[Tcstr1,Xcstr1]=ClassicalRungeKuttaFixedStep(@CSTR1, [t0, tf], T0, N, p);
%%
figs4 = tiledlayout(2,1);
figs4.Padding = 'compact';
figs4.TileSpacing = 'compact';
ax31 = nexttile;
plot(Tcstr1,Xcstr1)
pbaspect([5, 1, 1])
% xlim([0,2100])
grid on
ylabel('$$T[K]$$','interpreter','latex')
title("RK4 fixed step, "+"$$N=10^7$$",'interpreter','latex')
ax31 = nexttile;
plot(Tcstr1,Xcstr1)
pbaspect([5, 1, 1])
xlim([0,.001])
grid on
xlabel('$$t[s]$$','interpreter','latex')
ylabel('$$T[K]$$','interpreter','latex')
exportgraphics(figs4,'cstr1_RK4_n1e7.pdf','ContentType','vector')
% X in "XCV" is the factor multiplied to the Correct Volume used
disp('graph exported')

%% ADAPTIVE, CLassicalRungeKuttaAdaptiveStep.m
% Variables needed for adaptive step method:
abstol = 1e-9; % Usual / default value is 1e-6
reltol = 1e-9; % Default is 1e-3 in odeset.
% INITIAL STEP SIZE CALCULATION
%h0 = InitialStepSize(@CSTR1,t0,T0,p);
h0=0.0001;

[Tcstr1,Xcstr1]=ClassicalRungeKuttaAdaptiveStep(@CSTR1, [t0, tf], T0, h0, abstol, reltol, p);

fig4= figure(1);
plot(Tcstr1,Xcstr1)
% xlim([0 0.001])
xlabel('$$t[min]$$','interpreter','latex')
ylabel("Temperature")
grid on

%% Many flow rates ADAPTIVE, CLassicalRungeKuttaAdaptiveStep.m
% Variables needed for adaptive step method:
abstol = 1e-10; % Usual / default value is 1e-6
reltol = 1e-10; % Default is 1e-3 in odeset.
% INITIAL STEP SIZE CALCULATION
%h0 = InitialStepSize(@CSTR1,t0,T0,p);
h0=0.0001;

Ff(1) = 100;
Ff(2) = 250;
Ff(3) = 400;
Ff(4) = 550;
Ff(5) = 700;
Ff(6) = 850;
Ff(7) = 1000;
Tcstr=cell(1,7);
Xcstr=cell(1,7);

for k=1:7
    p.F=Ff(k);
    [Tcstr{k},Xcstr{k}]=ClassicalRungeKuttaAdaptiveStep(@CSTR1, [t0, tf], T0, h0, abstol, reltol, p);
end

%% s
figure('Position', [100 100 1000 500]);  % [left, bottom, width, height]
fig3A= figure(1);
for k=1:7
    plot(Tcstr{k}, Xcstr{k}, "DisplayName", "F="+string(Ff(k)) )
    hold on
end
hold off
pbaspect([5 1 1]);
lgd=legend("location", "eastoutside");
lgd.Title.String = '[mL/min]';

xlim([0,0.01])

ax=gca;
yticks = linspace(min(Xcstr{1}), max(Xcstr{1}), 4);
ax.YTick = yticks;
ax.YTickLabel = sprintfc('%.1e', yticks - 273.65);

title("CSTR1, RK4 Adaptive")
xlabel('$$t[min]$$','interpreter','latex')
ylabel('$$T[C]$$','interpreter','latex')
grid on
exportgraphics(fig3A,'cstr1_RK4_adapt_flows.pdf','ContentType','vector') % X in "XCV" is the factor multiplied to the Correct Volume used
disp('graph exported')