function xdot = CSTR3(t,x,p)
% CSTR 3 State 
%
%
xdot=zeros(3,1);
% Upnack
F = p.F;
V = p.V;
EaR=p.EaR;
k0=p.k0;
beta = p.beta;

cA=x(1);
cB=x(2);
T=x(3);

cAin = p.u(1);
cBin = p.u(2);
Tin  = p.u(3); % Boundary condition, in the 1 State case it happens to be the initial cond too.

kT= k0*exp(-EaR/T);
r = kT*cA*cB; %r=k(T)*cA*cB 
RA = -r;
RB = -2*r;
RT = beta*r;

cAdot = F/V * (cAin -cA) + RA;
cBdot = F/V * (cBin -cB) + RB;
cTdot = F/V *(Tin - T) + RT ;

xdot(1)=cAdot;
xdot(2)=cBdot;
xdot(3)=cTdot;