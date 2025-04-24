function Jac = VanDerPolJac(t,x,mu)
% JACVANDERPOL Jacobian for the Van der Pol Equation
%
% Syntax: Jac = JacVanDerPol(t,x,mu)
Jac = zeros(2,2); %(Row,column)
Jac(1,2) = 1.0;
Jac(2,1) = -2*mu*x(1)*x(2)-1.0;
Jac(2,2) = mu*(1-x(1)*x(1));
end