function Jac = CSTR3Jac(t,x,p)
% CSTR 3 State 
%

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
Tin  = p.u(3); % Boundary conditions, in the 1 State case this happens to be the initial cond too.

kT= k0*exp(-EaR/T);

% Jacobian
Jac = zeros(3,3);
Jac(1,1) = - F/V -kT*cA;
Jac(1,2) = -kT*cA;
Jac(1,3) = -kT*EaR*T^(-2);

Jac(2,1) = -2*kT*cB;
Jac(2,2) = -F/V - 2*kT*cA;
Jac(2,3) = -2*kT*EaR*T^(-2);

Jac(3,1) = beta*kT*cB;
Jac(3,2) = beta*kT*cA;
Jac(3,3) = - F/V +beta*kT*EaR*T^(-2);