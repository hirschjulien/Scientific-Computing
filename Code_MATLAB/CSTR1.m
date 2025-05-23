function xdot = CSTR1(t,x,p)
% CSTR 1 State where C=[T], therefore x, cAin, T are the same
%
%

% Upnack. This is unessecary (afaik) but it makes the code a bit clearer.
F = p.F;
V = p.V;
EaR=p.EaR;
k0=p.k0;
beta = p.beta;

cAin = p.u(1);
cBin = p.u(2);
Tin  = p.u(3); % Boundary condition, in the 1 State case it happens to be the initial cond too.

kT= k0*exp(-EaR/x);             % Arrheniious expression k(T)
cAT= cAin + (1/beta)*(Tin -x);  % Concentration A
cBT= cBin + (2/beta)*(Tin -x);  % Concentration B
r = kT*cAT*cBT;                 % Rate of reaction r
RT = beta*r;                    
Tdot = F/V *(Tin - x) + RT ;    % ODE, rate of change of temperature

xdot=Tdot;