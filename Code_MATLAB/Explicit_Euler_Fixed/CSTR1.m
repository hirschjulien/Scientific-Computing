function xdot = CSTR1(t,x,u,p)
% CSTR 1 State where C=[T], therefore x, cAin, T are the same
%
%

% Upnack
F = p.F;
V = p.V;
EaR=p.EaR;
k0=p.k0;
beta = p.beta;

cAin = u(1);
cBin = u(2);
Tin  = u(3); % Boundary condition, in the 1 State case it happens to be the initial cond too.

kT= k0*exp(-EaR/x);
cAT= cAin + (1/beta)*(Tin -x);
cBT= cBin + (2/beta)*(Tin -x);
r = kT*cAT*cBT;                 % HUSK!
RT = beta*r;
cAdot = F/V *(Tin - x) + RT ;

xdot=cAdot;