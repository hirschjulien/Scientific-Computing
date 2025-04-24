function [f,J] = PreyPredatorFunJac(t,x,a,b)
% The ode system
f = zeros(2,1);
f(1) = a*(1-x(2))*x(1);
f(2) = -b*(1-x(1))*x(2);
% The jacobian
J = zeros(2,2);
J(1,1) = a*(1-x(2));
J(2,1) = b*x(2);
J(1,2) = -a*x(1);
J(2,2) = -b*(1-x(1));
end