function [f, J] = VanDerPolFunJac(t,x,mu)

% JACVANDERPOL Jacobian for the Van der Pol Equation
%
% Syntax: J = JacVanDerPol(t,x,mu)
J = zeros(2,2);
j(1,1) = 0.0;
J(1,2) = 1.0;
J(2,1) = -2*mu*x(1)*x(2)-1.0;
J(2,2) = mu*(1-x(1)*x(1));

% VANDERPOL Implementation of the Van der Pol model
%
% Syntax: f = VanDerPol(t,x,mu)
f=zeros(2,1);
f(1) = x(2);
f(2) = mu*(1-x(1)*x(1))*x(2)-x(1);

end