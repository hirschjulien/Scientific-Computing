mu     = 3;
sigma  = 0.5;
x0     = [0.5; 0.5];
p      = [mu; sigma];

tf     = 5 * mu;
nw     = 1;
N      = 10000;
Ns     = 10;
seed   = 100;

[W, T, ~] = StdWienerProcess(tf, N, nw, Ns, seed);
X = zeros(length(x0), N+1, Ns);

for i = 1:Ns
    X(:,:,i) = ImplicitExplicit(...
        @VanderpolDrift, @VanderPolDiffusion1, ...
        T, x0, W(:,:,i), p);
end

% Deterministic solution, sigma = 0
Xd = ImplicitExplicit(@VanderpolDrift, @VanderPolDiffusion1, ...
    T, x0, W(:,:,i), [mu; 0.0]);

% --- Plotting ---
figure('Position', [100, 100, 1000, 500]);

% Plot x1(t)
subplot(2,1,1); hold on;
for i = 1:Ns
    plot(T, squeeze(X(1,:,i)), 'Color', [0.5, 0.5, 1]);
end
plot(T, Xd(1,:), 'k-', 'LineWidth', 2, 'DisplayName', 'Deterministic');
titleStr = sprintf(['ImplicitExplicit SDE Solver\n' ...
    '\\mu = %.1f, \\sigma = %.1f, x_0 = [%.1f, %.1f], Ns = %d'], ...
    mu, sigma, x0(1), x0(2), Ns);

title(titleStr, 'Interpreter', 'tex');
%xlabel('Time');
ylabel('x_1');

% Plot x2(t)
subplot(2,1,2); hold on;
for i = 1:Ns
    plot(T, squeeze(X(2,:,i)), 'Color', [1, 0.6, 0.6]);
end
plot(T, Xd(2,:), 'k-', 'LineWidth', 2, 'DisplayName', 'Deterministic');
%title('State x_2(t)');
xlabel('Time');
ylabel('x_2');

% Plot phase diagram
figure('Position', [100, 100, 600, 500]); hold on;
for i = 1:Ns
    plot(squeeze(X(1,:,i)), squeeze(X(2,:,i)), 'Color', [1, 0.6, 0.6]);
end
plot(Xd(1,:), Xd(2,:), 'k-', 'LineWidth', 2, 'DisplayName', 'Deterministic');
titleStr = sprintf(['ImplicitExplicit SDE Solver\n' ...
    '\\mu = %.1f, \\sigma = %.1f, x_0 = [%.1f, %.1f], Ns = %d'], ...
    mu, sigma, x0(1), x0(2), Ns);

title(titleStr, 'Interpreter', 'tex');
xlabel('x_1');
ylabel('x_2');


