function xdot = PFR1Fun(t,x,p)
%
%
%
% Unpacking variables
T   = x;
Tin = p.u(3);

n   = p.Nz;
dz  = p.dz;
nu  = p.nu;
% v  = p.v;          % Asks for velocity
F   = p.F;
A   = p.A;
v   = F/A;          % In our case "F" could be a vector  % HUSK!
DA  = p.DA;
k0  = p.k0;        
EaR = p.EaR;
beta = p.beta;
                
% Convection at finite volume interfaces
NconvA          = zeros(n+1,1);     %
NconvA(1,1)     = v*Tin;           % B.C. @ z=0
NconvA(2:n+1,1) = v*T(1:n,1);      % B.C. @ z=L

% Diffucsion at finite volume interfaces
JA              = zeros(n+1,1);     % Initialize data structure
JA(2:n,1)       = -DA/dz*(T(2:n,1) - T(1:n-1,1));     % ???

% Flux = convection + diffusion
NA              = NconvA + JA;

% Reaction and production rates in finte volumes
cAin    = p.u(1);
cBin    = p.u(2);
kT      = k0*exp(-EaR/T);         % HUSK! "T" is a vector
cAT     = cAin + (1/beta)*(Tin -T);
cBT     = cBin + (2/beta)*(Tin -T);
r       = kT*(cAT.*cBT);
RT      = beta*r;

% Differential equation (mass balances at finite volumes)
Tdot = (NA(2:n+1,1)-NA(1:n,1))/(-dz) + RT;             %
xdot = Tdot;                                           % Finalize output.