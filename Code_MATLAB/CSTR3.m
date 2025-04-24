function xdot = CSTR3(t,x,u,p)
% CSTR 1 State where C=[T], therefore x, cAin, T are the same
%
%

% Upnack
F = p.F;
V = p.V;
EaR=p.EaR;
k0=p.k0;
beta = p.beta;

cA=x(1);
cB=x(2);
cT=x(3);

cAin = u(1);
cBin = u(2);
Tin  = u(3); % Boundary condition, in the 1 State case it happens to be the initial cond too.

kT= k0*exp(-EaR/x);
r = kT*cA*cB; %r=k(T)*cA*cB 
RA = -r;
RB = -2*r;
RT = beta*r;

cAdot = F/V * (cAin -cA) + RA;
cBdot = F/V * (cBin -cB) + RB;
cTdot = F/V *(Tin - cT) + RT ;

xdot(1)=cAdot;
xdot(2)=cBdot;
xdot(3)=cTdot;