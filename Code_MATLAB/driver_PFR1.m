% driver for the PFR1Fun.m function
%
%
clear
clc
% --------------------
% Parameters
% --------------------
Nz = 100;                % [-]           Resolution of C in the spatial
% dimension
rho  = 1;              % [kg/m^3]        Density
c_p  = 4.186;            % [kj/(kg*K)]   Specific heat capacity
k0  = exp(24.6);        % [L/(mol*s)]   Arrhenius constant
EaR  = 8500.0;           % [K]           Activation energy
DHr  = -560.0;           % [kJ/mol]      Enthalpy of reaction
L    = 10;               % [m]           Length of reactor
A    = 0.1;              % [m^2]         Cross-sectional area of reactor
V    = A*L*1000;              % [m^3]         Reactor volume (=1m^3)
%D    = [0.1; 0.1; 0.1];  % [m^2/s]       Diffusion coefficients
D    = 0.1;
cAin = 1.6/2;            % [mol/L]       Concentration of A at Inlet
cBin = 2.4/2;            % [mol/L]       Concentration of B at Inlet
Tin  = 273.65;           % [K]           Temperature at Inlet (=0.5C)
% --------------------

% Depended parameters
beta = -DHr/(rho*c_p);
nu = [-1; -2; beta];
dz = L/Nz;

% Arbitrary vol. flow rate
F = 0.0083;              % Volumetric flowrate in 0.0083 [L/s] = 0.083 [m^3/s] = 500 mL/min, 0.5ml/min

% Create the parameter struct
p.nu    = nu;
p.k0    = k0;
p.dz    = dz;
p.Nz    = Nz;
p.EaR   = EaR;
p.DA    = D;
p.F     = F;
p.A     = A;
p.beta  = beta;
%v = F/A not p.v;
u = [cAin;cBin;Tin];    % Boundary conditions
p.u=u;
InitCond = zeros(Nz,1)+Tin;   % Initial conditions
t0 = 0;
tf = 60; % 35 minutes in seconds (=2100s)

% Set default font, font size for axes labels and titles
% set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultAxesTitleFontSizeMultiplier', 1.2);
set(groot, 'defaultAxesLabelFontSizeMultiplier', 1.2);

% (Optional) Set default line width and line style 
set(groot, 'defaultLineLineWidth', 1.8 );
%set(groot, 'defaultLineLineStyle', '--');

%% ODE 15s - No jacobian
clc
    
options = odeset('RelTol',1.0e-8,'AbsTol',1.0e-8);
[T,X]=ode15s(@PFR1Fun,[t0, tf],InitCond,options,p);

fig1= figure(2);
plot(T,X)
grid on

%%
x=linspace(1,100);
cmap=figure(1);
imagesc(x, T, X);
grid on
axis xy
cb = colorbar;         % create and get handle to colorbar
cb.Label.String = 'Temperature [K]';
hold on;
% Overlay contours
[x_grid, T_grid] = meshgrid(x, T);   % grid must match X size
contour(x_grid, T_grid, X, 10, 'k'); % 10 contour levels, black lines
hold off;
xlabel("Discrete space x")
ylabel('Time [s]')
title("PFR1 - oe15s, no Jacobian - tol:1e-8")

% exportgraphics(cmap,'pfr1_ode15s_noJac.pdf','ContentType','vector')
% disp('graph exported')


%% Implicit solver - ode15s
clc
options = odeset('Jacobian',@PFRJac,'RelTol',1.0e-6,'AbsTol',1.0e-6);
[T,X]=ode15s(@CSTR1, [t0, tf], InitCond, options, u, p);

fig2= figure(2);
plot (T,X)
grid on

%% RK4 fixed
clc
N=1e6;

[T,X]=ClassicalRungeKuttaFixedStep(@PFR1Fun, [t0, tf], InitCond, N, p);

%% RK4 fixed plot 
x=linspace(1,100);
cmap=figure(1);
imagesc(x, T, X);
grid on
axis xy
cb = colorbar;         % create and get handle to colorbar
cb.Label.String = 'Temperature [K]';
hold on;
% Overlay contours
[x_grid, T_grid] = meshgrid(x, T);   % grid must match X size
contour(x_grid, T_grid, X, 10, 'k'); % 10 contour levels, black lines
hold off;
xlabel("Discrete space x")
ylabel('Time [s]')
title("PFR1 - Rk4, fixed - tol:1e-8")

% exportgraphics(cmap,'pfr1_ode15s_noJac.pdf','ContentType','vector')
% disp('graph exported')