% CSTR1 Steady state T(F), C_A(T(F)), C_B(T(F))

clc
clear
% --------------------
% Parameters
% --------------------
p.rho  = 1.0;                   % [kg/L]        Density
p.c_p  = 4.186;                 % [kj/(kg*K)]   Specific heat capacity
p.DHr  = -560.0;                % [kJ/mol]      Enthalpy of reaction
p.beta = -p.DHr/(p.rho*p.c_p);
p.k0   = exp(24.6);             % [L/(mol*s)]   Arrhenius constant
p.EaR  = 8500.0;                % [K]           Activation energy
p.V    = 0.105;                 % [L]           For CSTR from Wahlgreen
p.DA   = [0.1; 0.1; 0.1];       % [m^2/s]       Diffusion coefficients
% p.nu = [-1; -2; p.beta];        % 3 state only!
% p.dz = p.L/p.Nz;                % PFR case only!
% Boundary conditions
cAin = 1.6/2;                   % [mol/L]       Concentration of A at Inlet
cBin = 2.4/2;                   % [mol/L]       Concentration of B at Inlet
Tin  = 273.65;                  % [C]           Temperature at Inlet (0.5C=273.65K)
u=[cAin;cBin;Tin];              % 3 state case
% --------------------
Abs0 = 273.15;
% the steady state temperature, Ts, chosen such that
% 0 ≤ CAs(Ts) ≤ CA,in and 0 ≤ CBs(Ts) ≤ CB,in
points = 100;
Ts  = linspace(0.5,100,points)+Abs0;
Fs  = zeros(1,points);
cAs = zeros(1,points);
cBs = zeros(1,points);

for j=1:points
    kT    = p.k0*exp(-p.EaR/Ts(j));
    cAs(j)   = cAin + (1/p.beta)*(Tin-Ts(j));
    cBs(j)   = cBin + (2/p.beta)*(Tin-Ts(j));
    r     = kT*cAs(j)*cBs(j);
    Fs(j) = p.V*p.beta*r/(Ts(j)-Tin);    
end

% Preparing variables for presentation
plotable_Fs = Fs*6e4; % multiplied by 6*10^4 to convert L/s to mL/min
plotable_Ts = Ts-Abs0;
%%
% Plotting
fig1=figure(1);
plot(plotable_Fs,plotable_Ts, 'LineWidth', 2)
pbaspect([4 1 1])
grid on
ylim([0,100]);
xlim([0,1000]);
xlabel('$$F_s [mL/min]$$','interpreter','latex', 'FontSize', 12)
ylabel('$$T_s[C]$$','interpreter','latex', 'FontSize', 12)
%exportgraphics(fig1,'cstr_Ts(Fs).pdf','ContentType','vector')
disp("Ολοκληρώθηκε!")

%%
clc
fig2=figure(2);
plot(plotable_Fs, cAs, "k-.", plotable_Fs, cBs, 'LineWidth', 2);
pbaspect([4 1 1])
legend({'$$C_As$$', '$$C_Bs$$'}, 'interpreter','latex', 'FontSize', 12, "Location", "southeast");
grid on
% ylim([0,100]);
xlim([0,1000]);
xlabel('$$F_s [mL/min]$$','interpreter','latex', 'FontSize', 12)
ylabel('$$C_s [mol/L]$$','interpreter','latex', 'FontSize', 12)
exportgraphics(fig2,'cstr_cs.pdf','ContentType','vector')
disp("Ολοκληρώθηκε!")