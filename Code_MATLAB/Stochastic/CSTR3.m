function [xdot,J] = CSTR3(t,x,u,p)
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

kT = k0 * exp(-EaR / cT);   % changed

r = kT*cA*cB; %r=k(T)*cA*cB 
RA = -r;
RB = -2*r;
RT = beta*r;

cAdot = F/V * (cAin -cA) + RA;
cBdot = F/V * (cBin -cB) + RB;
cTdot = F/V *(Tin - cT) + RT ;

xdot = [cAdot; cBdot; cTdot];

% Jacobian matrix (optional)
    if nargout > 1
        dkT_dT = kT * EaR / (cT^2);
        dr_dA = kT * cB;
        dr_dB = kT * cA;
        dr_dT = dkT_dT * cA * cB;

        dRA = -dr_dA;
        dRB = -2 * dr_dB;
        dRT = beta * dr_dT;

        J = zeros(3,3);
        J(1,1) = -F/V + dRA;
        J(1,2) = 0 + (-kT * cA);
        J(1,3) = 0 + (-dkT_dT * cA * cB);

        J(2,1) = 0 + (-2 * kT * cB);
        J(2,2) = -F/V + dRB;
        J(2,3) = 0 + (-2 * dkT_dT * cA * cB);

        J(3,1) = 0;
        J(3,2) = 0;
        J(3,3) = -F/V + dRT;
    end

%xdot(1)=cAdot;
%xdot(2)=cBdot;
%xdot(3)=cTdot;