function xdot = PFR3Fun(t,x,p)
%
%
%
% Unpacking variables
x
cA = x(:,1);
cB = x(:,2);
T  = x(:,3);
cAin = p.u(1);
cBin = p.u(2);
Tin  = p.u(3);

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
kT      = k0*exp(-EaR/T);         % HUSK! "T" is a vector
r       = kT*(cA.*cB);
R      = nu*r;

% Differential equation (mass balances at finite volumes)
Tdot(1) = (NA(2:n+1,1)-NA(1:n,1))/(-dz) + R(1);    
Tdot(2) = (NA(2:n+1,1)-NA(1:n,1))/(-dz) + R(2);    
Tdot(3) = (NA(2:n+1,1)-NA(1:n,1))/(-dz) + R(3);             %
xdot = Tdot;                                           % Finalize output.