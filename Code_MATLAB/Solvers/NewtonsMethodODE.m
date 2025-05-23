function x = NewtonsMethodODE(FunJac, tk, xk, dt, xinit, tol, maxit, varargin)
% Inputs:
    % FunJac : a function that returns both the function to be solved and
        % its jacobian
    % tk := the time at the specific time-step we are solving
    % xk := the value at such time
    % dt := the time-step
    % xinit := an initial guess of the value we are looking for
    % tol := tolerance, 
    % maxit := maximum number of iterations. These laste two steps resemble
        % memory allocation and are here to keep the computational cost under
        % control
    % varargin : rest arguments and function parametres

k=0;
t = tk + dt;
x = xinit;
[f,J] = feval(FunJac, t, x, varargin{:});
R = x - f*dt - xk;      % The residual function
I = eye(length(xk));
while( (k<maxit) && (norm(R, 'inf')>tol) )
    k = k+1;            % Can also be used to count number of function calls
    dRdx = I - J*dt;        % We also set this equaiton as "M" in the slides
    dx = dRdx\R;
    x = x - dx;
    [f, J] = feval(FunJac, t, x, varargin{:});
    R = x - f*dt - xk;
end

