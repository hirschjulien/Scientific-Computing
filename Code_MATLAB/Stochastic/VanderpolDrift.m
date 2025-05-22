function [f, J] = VanderpolDrift(t, x, p)
    mu = p(1);
    tmp = mu * (1.0 - x(1) * x(1));
    f = zeros(2,1);
    f(1) = x(2);
    f(2) = tmp * x(2) - x(1);

    if nargout > 1
        J = [0 1; -2 * mu * x(1) * x(2) - 1.0, tmp];
    end
end
