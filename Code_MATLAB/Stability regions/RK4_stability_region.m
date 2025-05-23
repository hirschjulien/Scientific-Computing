% Stability region for RK4


[X, Y] = meshgrid(-5:.01:1, -3:.01:3); 
z = X + 1i*Y; 
R = 1 + z + .5*z.^2 + (1/6)*z.^3 + (1/24)*z.^4; 
Rhat = abs(R); 
Rhat = Rhat.*(Rhat<1);  %# here truncation occurs

fig1 = figure(1);
imagesc([min(X(:)) max(X(:))],[min(Y(:)) max(Y(:))], Rhat) 
colormap(flipud(bone))
cb = colorbar();
cb.Title.String = "R(z)";
xlabel('Real(z)')
ylabel('Imag(z)')
grid on;
exportgraphics(fig1,'rk4_stability_region.pdf','ContentType','vector')
disp('fin')

%Reference: https://www.physicsforums.com/threads/plotting-runge-kutta-4-stability-region.1023646/

