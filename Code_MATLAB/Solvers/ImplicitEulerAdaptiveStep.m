function [T,X] = ImplicitEulerAdaptiveStep(FunJac, tspan, x0, h0, abstol, reltol, varargin)

% h_{k+1} = [ (ε/r_{k+1})^(1/2) ]*h_k

% Error control parameters
epstol = 0.8;   % Safety factor.  Also shown as α in slides. We choose 0.8 instead of 1 in case we have inaccuracies at ch^2 calculation of |e|. A simple controller of the numerical process. 
facmin = 0.1;   % Maximum decrease factor. Heuristic
facmax = 5.0;   % Maximum increase factor. Heuristic
% Integration interval
t0 = tspan(1);  % Initial time
tf = tspan(2);  % Final time
% Initial conditions
t = t0;
h = h0;
x = x0;
% Output
T = t;
X = x';
% ImplicitEuler parameters
tol = 1.0e-8;
maxit = 100; 
%% Main algorithm

while t<tf
    if t+h>tf
        h = tf-t;
    end
    f=feval(FunJac,t,x,varargin{:});

    AcceptStep = false;
    while ~AcceptStep
        % Take step of size h
        xinit = x + h*f; % Our initial guess of the function value at this time=t_k
        x1 = NewtonsMethodODE(FunJac, t, x, h, xinit, tol, maxit, varargin{:});

        % Take step of size h/2
        hm = 0.5*h;
        tm = t + hm;
        %xm = x + hm*f;  % ANOTHER EULER STEP
        xinitm = x + hm*f; % Our initial guess of the function value at this time=t_k
        xm = NewtonsMethodODE(FunJac, t, x, hm, xinitm, tol, maxit, varargin{:});
        fm = feval(FunJac, tm, xm, varargin{:}); % Function evaluated with half step
        x1hat = NewtonsMethodODE(FunJac, tm, xm, hm, xinit, tol, maxit, varargin{:});

        % Error estimation
        e = x1hat - x1;
        r = max( abs(e)./max( abstol, abs(x1hat).*reltol ) );   
            % For all components in the vector, we can consider it the 
            % ρ(e_{n+1}) function

        AcceptStep = (r <= 1.0);    % Checks if the accuracy of the step is acceptable
        if AcceptStep
            t = t+h;        % If the step is not accepted, the initial t and x are kept
            x = x1hat;

            T = [T;t];
            X = [X;x'];
        end
        % Asymptotic step size controller
        h = max(facmin, min(sqrt(epstol/r), facmax))*h;
            % If the error is rly rly small we should increase the step
            % size but we also control the increase and deacrease.
    end
end




