function Jac = CSTR1Jac(t,x,p)
% CSTR 1 State where C=[T], therefore x, cAin, T are the same
%
%

% Upnack
F = p.F;
V = p.V;
EaR=p.EaR;
k0=p.k0;
beta = p.beta;

cAin = p.u(1);
cBin = p.u(2);
Tin  = p.u(3); % Boundary condition, in the 1 State case it happens to be the initial cond too.

kT= k0*exp(-EaR/x);
cAT= cAin + (1/beta)*(Tin -x);
cBT= cBin + (2/beta)*(Tin -x);
Jac = -F/V +kT*(-cBT -2*cAT + beta*EaR*x^(-2));