% Define complex plane grid
th = linspace(-5, 5, 400);
[X, Y] = meshgrid(th);
Z = X + 1i*Y;

% Compute magnitude of the amplification factor
R = 1 ./ (1 - Z);
abs_R = abs(R);

% Plot region where |R(z)| < 1
figure
contourf(X, Y, abs_R, [0 1], 'LineColor','none')  % Stability region
%colormap([0.5 0.8 1])  % Light blue
hold on
%contour(X, Y, abs_R, [1 1], 'k', 'LineWidth', 2)    % Boundary
xlabel('Real(z)')
ylabel('Imag(z)')
title('Absolute Stability Region of Implicit Euler Method')
grid on
axis equal

%%
% Define complex plane grid
x = linspace(-5, 5, 400);
y = linspace(-5, 5, 400);
[X, Y] = meshgrid(x, y);
Z = X + 1i*Y;

% Compute magnitude of the amplification factor
R = 1 ./ (1 - Z);
abs_R = abs(R);
% Logical mask for stable region
stable_region = abs_R < 1;

% Plot
figure
hold on

% Plot the stable region (light blue), rest is white by default
contourf(X, Y, stable_region, [1 1], 'FaceColor', [0.7 0.85 1], 'LineColor', 'none')

% Plot boundary (|R(z)| = 1)
contour(X, Y, abs_R, [1 1], 'k', 'LineWidth', 1)

% Axis styling
xlabel('Real(z)')
ylabel('Imag(z)')
title('Absolute Stability Region of Implicit Euler Method')
axis equal
axis([-5 5 -5 5])  % optional: fix axis range
grid on
set(gca, 'Layer', 'top')      % Keep grid/ticks above patch
box on

%%
[X, Y] = meshgrid(-3:.01:3, -3:.01:3); 
z = X + 1i*Y; 
R = 1./(1-z); 
Rhat = abs(R); 
Rhat = Rhat.*(Rhat<1);  %# here I truncate 

fig1 = figure(1);
imagesc([min(X(:)) max(X(:))],[min(Y(:)) max(Y(:))], Rhat) 
colormap(flipud(bone))
cb = colorbar();
cb.Title.String = "R(z)";
xlabel('Real(z)')
ylabel('Imag(z)')
grid on;
exportgraphics(fig1,'implicit_euler_stability_region.pdf','ContentType','vector')
disp('fin')
%Reference: https://www.physicsforums.com/threads/plotting-runge-kutta-4-stability-region.1023646/