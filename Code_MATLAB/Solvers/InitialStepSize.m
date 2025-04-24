function h0 = InitialStepSize(fun,t,x,varargin)
d0 = norm(x);
d1 = norm(feval(fun,t,x,varargin{:}));
if any([d0, d1] < 1e-5)
    h0 = 1e-6;
else
    h0 = .01*(d0/d1);
end