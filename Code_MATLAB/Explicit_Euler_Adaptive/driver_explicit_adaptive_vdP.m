%%% Predator Prey
%tspan = [0,50];
%p = [1,1];
%x0 = [2;2];
%h0 = 1;

%%% van der Pol
tspan = [0,50];
mu = 10;
x0 = [2;2];
h0 = 1;

reltol = [1e-1, 1e-3, 1e-5, 1e-8];
abstol = 1e-5;

figure('Position', [100, 100, 1000, 600]);

% Prey subplot
subplot(2,1,1);
hold on;  % Keep all plots in same axes
for tol = reltol
    [T,X] = ExplicitEulerAdaptiveStep(@VanDerPol, tspan, x0, h0, abstol, tol, mu);
    plot(T, X(:,1), 'DisplayName', sprintf('reltol = %d', tol));
end
[Tcorrect,Xcorrect] = ode15s(@(t,x) VanDerPol(t,x,mu), tspan, x0);
plot(Tcorrect, Xcorrect(:,1), 'DisplayName', sprintf('ode15s'));

xlabel('Time');
ylabel('X(:,1)');
legend;
hold off;

% Predator subplot
subplot(2,1,2);
hold on;
for tol = reltol
    [T,X] = ExplicitEulerAdaptiveStep(@VanDerPol, tspan, x0, h0, abstol, tol, mu);
    plot(T, X(:,2), 'DisplayName', sprintf('reltol = %d', tol));
end
[Tcorrect,Xcorrect] = ode15s(@(t,x) VanDerPol(t,x,mu), tspan, x0);
plot(Tcorrect, Xcorrect(:,2), 'DisplayName', sprintf('ode15s'));

xlabel('Time');
ylabel('X(:,2)');
legend;
hold off;
