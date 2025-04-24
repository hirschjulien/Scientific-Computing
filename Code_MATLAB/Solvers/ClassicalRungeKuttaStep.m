function [t1,x1] = ClassicalRungeKuttaStep(fun, t, x, f, h, varargin)

h2 = 0.5*h;
alpha = h/6; % from the RK Buthcer tableu
beta = h/3;

T1 = t;
X1 = x;
F1 = f;

T2 = x + h2;
X2 = x + h2*F1;
F2 = feval(fun, T2, X2, varargin{:});

T3 = T2;
X3 = x + h2*F2;
F3 = feval(fun, T3, X3, varargin{:});

T4 = t+h;
X4 = x + h*F3;
F4 = feval(fun, T3, X4, varargin{:});

t1 = T4;
x1 = x + alpha*(F1+F4) + beta*(F2+F3); % Again, from the Classical RK Butcher tableu

end