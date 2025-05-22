function g = PredatorPreyDiffusion_Dependent(t, x, p)
    sigma = p(3);
    g = [sigma * (1 + x(1)^2); sigma * (1 + x(2)^2)];
end
