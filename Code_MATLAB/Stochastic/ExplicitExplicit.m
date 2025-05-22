function X = ExplicitExplicit(ffun, gfun, T, x0, W, varargin)

N  = length(T);
nx = length(x0);
X  = zeros(nx, N);

X(:,1) = x0;
for k = 1:N-1
    f  = feval(ffun, T(k), X(:,k), varargin{:});
    g  = feval(gfun, T(k), X(:,k), varargin{:});
    dt = T(k+1) - T(k);
    dW = W(:,k+1) - W(:,k);
    psi = X(:,k) + g * dW;
    X(:,k+1) = psi + f * dt;
end

end
