% driver_CSTR3 Script
%
%
%clc
%clear
% --------------------
% Parameters
% --------------------
clear p  % Remove any conflicting definition
p = struct();  % Explicitly define p as a struct
% p.Nz = 100;                     % [-]           Resolution of C in the spatial dimension
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
p.nu = [-1; -2; p.beta];        % 3 state only!
% p.dz = p.L/p.Nz;                % PFR case only!
p.sigma = 0.1;
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
InitCond = [0;0;Tin];       % Initial conditions

t0 = 0;
tf = 5;      % Data point every 10 seconds
nw     = 1;
N      = 10000;
Ns     = 10;
seed   = 100;

[W, T, ~] = StdWienerProcess(tf, N, nw, Ns, seed);

% aa = linspace(0, 1000, 2500);  % X-axis
% mu = 500;                     % Mean
% sigma = 100;                  % Standard deviation
% p.F = (1/(sigma * sqrt(2*pi))) * exp(-0.5 * ((aa - mu)/sigma).^2);

ffun = @(t, x, p, u) CSTR3(t, x, u, p);  % existing function
gfun = @(t, x, p, u) [0; 0; (p.F / p.V) * p.sigma];  % noise only in temperature

for i = 1:Ns
    W_i = squeeze(W(:,:,i));  % size: [1, N+1]
    X(:,:,i) = ImplicitExplicit(ffun, gfun, T, InitCond, W_i, p, u);
end

% Solution without stochastic influence i.e. sigma=0
p_d = p;
p_d.sigma = 0;

[W_dummy, ~, ~] = StdWienerProcess(tf, N, 1, 1, 999);  % dummy Wiener input
W0 = squeeze(W_dummy(:,:,1));
Xd = ImplicitExplicit(ffun, gfun, T, InitCond, W0, p_d, u);  % same functions, zero sigma


figure;
titles = {'C_A', 'C_B', 'T'};

colors = lines(Ns);  % Ns distinct colors

for j = 1:3  % for each state variable
    subplot(3,1,j); hold on;
  
    for i = 1:Ns
        plot(T, squeeze(X(j,:,i)), 'Color', colors(i,:), 'LineWidth', 1);
    end
    
 % Plot deterministic reference
    plot(T, Xd(j,:), 'k-', 'LineWidth', 2, 'DisplayName', 'Deterministic');

    ylabel(titles{j});
end

xlabel('Time [min]');
sgtitle('Stochastic CSTR3 with temperature noise only, implicit-explicit');

