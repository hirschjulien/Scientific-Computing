function g = PredatorPreyDiffusion1(t, x, p)
    sigma = p(3);  % same sigma for both
    g = [sigma; sigma];
end
