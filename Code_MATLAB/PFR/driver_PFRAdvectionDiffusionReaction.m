% driver for the PFRAdvectionDiffusionReaction.m functionm
%
%
clear
clc
% --------------------
% Parameters
% --------------------
Nz = 10;                % [-]           Resolution of C in the spatial
% dimension
rho  = 1.0;              % [kg/L]        Density
c_p  = 4.186;            % [kj/(kg*K)]   Specific heat capacity
k0  = exp(24.6);        % [L/(mol*s)]   Arrhenius constant
EaR  = 8500.0;           % [K]           Activation energy
DHr  = -560.0;           % [kJ/mol]      Enthalpy of reaction
L    = 10;               % [m]           Length of reactor
A    = 0.1;              % [m^2]         Cross-sectional area of reactor
V    = A*L;              % [m^3]         Reactor volume (=1m^3)
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
F = 0.0166;              % Volumetric flowrate in [m^3/s] 83e-3 = 500 mL/min

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
InitCond = zeros(Nz,1)+Tin;   % Initial conditions

% ODE 4-5
clc
t0 = 0;
tf = 35;     % 35 minutes in seconds (=2100s)
options = odeset('RelTol',1.0e-3,'AbsTol',1.0e-3);
[T,X]=ode45(@PFRAdvectionDiffusionReaction,[t0, tf],InitCond,options,u(3),p);


fig1= figure(1);
plot(T,X)
grid on

%% Implicit solver - ode15s
