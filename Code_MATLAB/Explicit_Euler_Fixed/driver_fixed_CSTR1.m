% driver_CSTR1 Script
%
%
%clc
%clear
% --------------------
% Parameters
% --------------------
p.Nz = 100;                     % [-]           Resolution of C in the spatial dimension
p.rho  = 1.0;                   % [kg/L]        Density
p.c_p  = 4.186;                 % [kj/(kg*K)]   Specific heat capacity
p.k0   = exp(24.6);             % [L/(mol*s)]   Arrhenius constant
p.EaR  = 8500.0;                % [K]           Activation energy
p.DHr  = -560.0;                % [kJ/mol]      Enthalpy of reaction
p.L    = 10;                    % [m]           Length of reactor
p.A    = 0.1;                   % [m^2]         Cross-sectional area of reactor
p.V    = p.A*p.L*1000;          % [L]           Reactor volume (A*L=1m^3)(1m^3=1000L)
p.DA   = [0.1; 0.1; 0.1];       % [m^2/s]       Diffusion coefficients
p.beta = -p.DHr/(p.rho*p.c_p);
% p.nu = [-1; -2; p.beta];        % 3 state only!
% p.dz = p.L/p.Nz;                % PFR case only!
% Boundary conditions
cAin = 1.6/2;                   % [mol/L]       Concentration of A at Inlet
cBin = 2.4/2;                   % [mol/L]       Concentration of B at Inlet
Tin  = 273.65;                  % [K]           Temperature at Inlet (=0.5C)
u=[cAin;cBin;Tin];            % 3 state case
% --------------------
 p.F=650;                        % [mL/min]     Volumetric Flowrate

% Create the parameter struct
%v = F/A not p.v;
% u = Tin;      % Boundary condition
T0 = Tin;       % Initial conditions

t0 = 0;
tf = 35*6;      % Data point every 10 seconds

% aa = linspace(0, 1000, 2500);  % X-axis
% mu = 500;                     % Mean
% sigma = 100;                  % Standard deviation
% p.F = (1/(sigma * sqrt(2*pi))) * exp(-0.5 * ((aa - mu)/sigma).^2);

tspan = [0, 5];
N_values = [10, 50, 100, 500];

% --- Plotting ---
figure('Position', [100, 100, 800, 400]);
hold on;

for N = N_values
    [T, X] = ExplicitEulerFixedStepSize(@CSTR1, tspan, N, Tin, u, p);
    plot(T, X, 'DisplayName', sprintf('N = %d', N));
end

% Reference solution using ode15s
[Tref, Xref] = ode45(@(t,x) CSTR1(t,x,u,p), tspan, Tin);
plot(Tref, Xref, 'k--', 'DisplayName', 'ode45');

xlabel('Time [min]');
ylabel('Temperature T(t) [K]');
title('CSTR1: Temperature evolution');
legend;
hold off;