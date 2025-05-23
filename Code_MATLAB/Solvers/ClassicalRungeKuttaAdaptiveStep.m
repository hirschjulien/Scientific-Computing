function [T,X] = ClassicalRungeKuttaAdaptiveStep( ...
    fun, tspan, x0, h0, abstol, reltol, varargin)
% varargin{:} expands the contents of varargin so that each element is
% passed as a separate input

% Error controller
epstol = 0.8;   % Target accuracy. We choose 0.8 instead of 1 in case we have inaccuracies at ch^2 calculation of |e|. A simple controller of the numerical process. 
kpow = 0.2;
facmin = 0.1;   % Maximum decrease factor
facmax = 5.0;   % Maximum increase factor
% Integration interval
t0 = tspan(1);
tf = tspan(2);
% Initial Conditions
t = t0;
h = h0;
x = x0;
% Output
T = t;
X = x';

%% Main Algorithm
while t < tf
    if (t+h > tf)
        h = tf - t;
    end
    f = feval(fun, t, x, varargin{:});

    AcceptStep = false;
    while ~AcceptStep       % "~" is MATLAB's NOT logival operator
        % Take step of size h
        [t1,x1] = ClassicalRungeKuttaStep(fun,t,x,f,h,varargin{:});

        % Take step of size h/2
        hm = 0.5*h;         % Halfing the step size
        [tm,xm]        = ClassicalRungeKuttaStep(fun, t,x,f,hm,varargin{:});
        fm             = feval(fun,tm,xm,varargin{:});
        [t1hat, x1hat] = ClassicalRungeKuttaStep(fun,tm,xm,fm,hm,varargin{:});
            % Here t1,x1 are re-calculated by a double step and below
            % accuracy is examined
        
        % Error estimation
        e = x1hat-x1;
        r = max( abs(e)./max(abstol, abs(x1hat).*reltol) );

        AcceptStep = (r<=1.0);
        if AcceptStep
            % If the double-step calculation is found more accurate, its
            % results are used.
            t=t+h;
            x=x1hat;
            T=[T;t];
            X=[X;x'];
        end 
        % Assymptotic step size controller
        h = max( facmin, min( (epstol/r)^kpow, facmax ) )*h;
    end     % End of nested while loop
end         % End of initial while loop







