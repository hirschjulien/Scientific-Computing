function [f, J] = PredatorPreyDrift(t, x, p)
    a = p(1);
    b = p(2);

    x1 = x(1);
    x2 = x(2);

    % Drift vector
    f = zeros(2,1);
    f(1) = a * (1 - x2) * x1;
    f(2) = -b * (1 - x1) * x2;

    % Jacobian matrix
    if nargout > 1
        J = [ a*(1 - x2),     -a*x1;
             -b*x2,           b*(1 - x1)];
    end
end
