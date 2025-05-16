function xdot = PFR3Fun(t,x,u,p)
%
%
%
% Unpacking variables
cA = x(:,1);
cB = x(:,2);
T  = x(:,s3);
cAin = u(1);
cBin = u(2);
Tin  = u(3);

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
NconvA          = zeros(n+1,3);    %
NconvA(1,1)     = v*cAin;           % B.C. @ z=0
NconvA(1,2)     = v*cBin;           % B.C. @ z=0
NconvA(1,3)     = v* Tin;           % B.C. @ z=0
NconvA(2:n+1,1) = v*cA(1:n,1);      % B.C. @ z=L
NconvA(2:n+1,2) = v*cB(1:n,1);      % B.C. @ z=L
NconvA(2:n+1,3) = v* T(1:n,1);      % B.C. @ z=L

% Diffucsion at finite volume interfaces
JA              = zeros(n+1,3);     % Initialize data structure
JA(2:n,1)       = -DA/dz*(cA(2:n,1) - cA(1:n-1,1));     % ???
JA(2:n,2)       = -DA/dz*(cB(2:n,2) - cB(1:n-1,2));     % ???
JA(2:n,3)       = -DA/dz*( T(2:n,3)  - T(1:n-1,3));     % ???

% Flux = convection + diffusion
NA              = NconvA + JA;

% Reaction and production rates in finte volumes
cAin    = u(1);
cBin    = u(2);
kT      = k0*exp(-EaR/T);         % HUSK! "T" is a vector
r       = kT*(cA.*cB);
RT      = beta*r;

% Differential equation (mass balances at finite volumes)
Tdot = (NA(2:n+1,1)-NA(1:n,1))/(-dz) + RA;    
Tdot = (NA(2:n+1,1)-NA(1:n,1))/(-dz) + RB;    
Tdot = (NA(2:n+1,1)-NA(1:n,1))/(-dz) + RT;             %
xdot = Tdot;                                           % Finalize output.