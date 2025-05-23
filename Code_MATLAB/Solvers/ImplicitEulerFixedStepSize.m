    function [T,X] = ImplicitEulerFixedStepSize(FunJac,tspan,N,tol,x0,p)
% uhmm

% Unpack
t0 = tspan(1);
tN = tspan(2);

% "Compute step size and allocate memory"
dt = (tN-t0)/N;
nx = size(x0,1);
X = zeros(nx,N+1);
T = zeros(1,N+1); 

%tol = 1.0e-8;
maxit = 100; 

% Eulers implicit method
T(:,1) = t0; % Initialise time
X(:,1) = x0; % Initialise function value
for k=1:N
    f = feval(FunJac, T(k), X(:,k), p);
    T(:,k+1) = T(:,k) + dt;
    xinit = X(:,k) + f*dt; % Our initial guess of the function value at this time=t_k
    X(:,k+1) = NewtonsMethodODE(FunJac, T(:,k), X(:,k), dt, xinit, tol, maxit, p);
end

% For MATLAB compatibility
T=T';
X=X';

end