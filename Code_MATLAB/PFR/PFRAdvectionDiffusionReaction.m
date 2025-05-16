function xdot = PFRAdvectionDiffusionReaction(t,x,u,p)
%
%
%
cA   = x;
cAin = u;

n   = p.Nz;
dz  = p.dz;

% v  = p.v;          % Asks for velocity
v   = p.F/p.A;          % In our case it depends on "F" which is a vector  % HUSK!
DA  = p.DA;
% k = p.k           % Asks for k(T)
k0  = p.k0;        
EaR = p.EaR;
k   = k0*exp(-EaR/cA);                   % HUSK! "T" is a vector
beta= p.beta;

% Convection at finite volume interfaces
NconvA          = zeros(n+1,1);     %
NconvA(1,1)     = v*cAin;           % B.C. @ z=0
NconvA(2:n+1,1) = v*cA(1:n,1);      % B.C. @ z=L

% Diffucsion at finite volume interfaces
JA              = zeros(n+1,1);     % Initialize data structure
JA(2:n,1)       = (-DA/dz)*(cA(2:n,1) - cA(1:n-1,1));     % ???

% Flux = convection + diffusion
NA = NconvA + JA;

% Reaction and production rates in finte volumes
r = k*(cA); 
RA = beta*r;
% % OR
% r = k*cA;
% RA = -nu*r;
% % OR
% r = k*(cA.*cA);                   % For r=k(T)C_a^2


% Differential equation (mass balances at finite volumes)
cAdot = (NA(2:n+1,1)-NA(1:n,1))/(-dz) + RA;             %
xdot = cAdot;                                           % Finalize output.