function [T,X] = ClassicalRungeKuttaFixedStep(fun, tspan, x0, N, varargin)
% Requires file:    ClassicalRungeKuttaStep.m
% varargin{:} expands the contents of varargin so that each element is
% passed as a separate input

% Integration interval
t0 = tspan(1);
tf = tspan(2);
% Initial Conditions
t = t0;
h=(tf-t0)/N;
x = x0;
% Output
T(:,1) = t0;
X(:,1) = x0;

for k=1:N
    f = feval(fun, T(k), X(:,k), varargin{:});
    [~,x1] = ClassicalRungeKuttaStep(fun, T(k), X(:,k), f, h, varargin{:});
    T(:,k+1) = T(:,k) + h;
    X(:,k+1) = x1;
end

T=T';
X=X';
end

