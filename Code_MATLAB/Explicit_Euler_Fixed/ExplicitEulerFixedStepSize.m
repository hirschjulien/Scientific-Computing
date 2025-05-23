function [T,X] = ExplicitEulerFixedStepSize(fun,tspan,N,x0,varargin)

ta = tspan(1);
tb = tspan(2);
% Compute step size and allocate memory
dt = (tb-ta)/N;
nx = size(x0,1);
X = zeros(nx,N+1);
T = zeros(1,N+1);

% Eulers Explicit Method
T(:,1) = ta;
X(:,1) = x0;
for k=1:N
    f = feval(fun,T(k),X(:,k),varargin{:});
    T(:,k+1) = T(:,k) + dt;
    X(:,k+1) = X(:,k) + f*dt;
end

% Form a nice table for the result
T = T';
X = X';
