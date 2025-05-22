sigma  = 0.05;
x0     = [2; 2];
p      = [1; 1; sigma];

tf     = 20;
nw     = 1;
N      = 10000;
Ns     = 10;
seed   = 100;

[W, T, ~] = StdWienerProcess(tf, N, nw, Ns, seed);
X = zeros(length(x0), N+1, Ns);

for i = 1:Ns
    X(:,:,i) = ImplicitExplicit(...
        @PredatorPreyDrift, @PredatorPreyDiffusion2, ...
        T, x0, W(:,:,i), p);
end

% Deterministic solution, sigma = 0
Xd = ImplicitExplicit(@PredatorPreyDrift, @PredatorPreyDiffusion2, ...
    T, x0, W(:,:,i), [1; 1; 0.0]);

% --- Plotting ---
figure('Position', [100, 100, 1000, 500]);

% Plot x1(t)
subplot(2,1,1); hold on;
for i = 1:Ns
    plot(T, squeeze(X(1,:,i)), 'Color', [0.5, 0.5, 1]);
end
plot(T, Xd(1,:), 'k-', 'LineWidth', 2, 'DisplayName', 'Deterministic');
titleStr = sprintf(['ImplicitExplicit SDE Solver\n' ...
    '\\sigma = %.2f, x_0 = [%.1f, %.1f], Ns = %d'], ...
    sigma, x0(1), x0(2), Ns);

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
titleStr = sprintf(['ExplicitExplicit SDE Solver\n' ...
    '\\sigma = %.2f, x_0 = [%.1f, %.1f], Ns = %d'], ...
    sigma, x0(1), x0(2), Ns);

title(titleStr, 'Interpreter', 'tex');
xlabel('x_1');
ylabel('x_2');


