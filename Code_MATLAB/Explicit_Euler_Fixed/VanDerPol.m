function xdot = VanDerPol(t,x,mu)
% VANDERPOL Implementation of the Van der Pol model
%
% Syntax: xdot = VanDerPol(t,x,mu)
xdot=zeros(2,1);
xdot(1) = x(2);
xdot(2) = mu*(1-x(1)*x(1))*x(2)-x(1);
end




