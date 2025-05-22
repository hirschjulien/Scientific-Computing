%%% Predator Prey
%tspan = [0,50];
%N_values = [1000, 5000, 10000, 50000];
%p = [1,1];
%x0 = [2;2];

%%% Van der Pol
tspan = [0,50];
N_values = [1000, 5000, 10000, 50000];
mu = 3;
x0 = [2;0];

figure('Position', [100, 100, 1000, 600]);

% Prey subplot
subplot(2,1,1);
hold on;  % Keep all plots in same axes
for n = N_values
    [T,X] = ExplicitEulerFixedStepSize(@VanDerPol, tspan, n, x0, mu);
    plot(T, X(:,1), 'DisplayName', sprintf('N = %d', n));
end
[Tcorrect,Xcorrect] = ode15s(@(t,x) VanDerPol(t,x,mu), tspan, x0);
plot(Tcorrect, Xcorrect(:,1), 'DisplayName', sprintf('ode45'));

xlabel('Time');
ylabel('X(:,1)');
legend;
hold off;

% Predator subplot
subplot(2,1,2);
hold on;
for n = N_values
    [T,X] = ExplicitEulerFixedStepSize(@VanDerPol, tspan, n, x0, mu);
    plot(T, X(:,2), 'DisplayName', sprintf('N = %d', n));
end
[Tcorrect,Xcorrect] = ode15s(@(t,x) VanDerPol(t,x,mu), tspan, x0);
plot(Tcorrect, Xcorrect(:,2), 'DisplayName', sprintf('ode45'));

xlabel('Time');
ylabel('X(:,2)');
legend;
hold off;
