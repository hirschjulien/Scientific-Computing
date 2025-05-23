%% Prey Predator explicit solver (non-stiff system)
% Also known as the Lotka-Volterra model in biotechnology
clear
clc

a = 1;
b = 1;
x0 = [2; 2];
options = odeset('RelTol',1.0e-6,'AbsTol',1.0e-6);
[T,X]=ode45(@PreyPredator,[0 50],x0,options,a,b);

plotter_prey_pred(T,X,a,b,"ode45")

% % Plots
% fig1 = figure(1);
% plot(T,X(:,1));
% title('')
% xlabel('time') 
% ylabel('x_1(t)') 
% fig2 = figure(2);
% plot(T,X(:,2));
% xlabel('time') 
% ylabel('x_2(t)') 
% fig3 = figure(3);
% plot(X(:,1),X(:,2));
% title('Prey - Predator population cycle')
% xlabel('x_1(t)') 
% ylabel('x_2(t)') 
%% implicit solver (implementation as IF it were a STIFF system)
% Requires the PrePredator.m function as well as its Jacobian
% JacPreyPredator.m function.
clear
clc

a = 1;  % Model
b = 42;  % parametres
x0 = [2; 2]; % Initial values, here the population of prey and predator
options = odeset('Jacobian',@PreyPredatorJac,'RelTol',1.0e-6,'AbsTol',1.0e-6); % Now the Jacobian is introduced in the solver options
[T,X]=ode15s(@PreyPredator,[0 50],x0,options,a,b); % The solver is changed to ode15s

plotter_prey_pred(T,X,a,b,"ode15s")

% figs2 = tiledlayout(2,2);
% figs2.Padding = 'compact';
% figs2.TileSpacing = 'compact';
% nexttile
% plot(T,X(:,1));
% title('One population over time');
% xlabel('time') 
% ylabel('x_1(t)') 
% nexttile
% plot(T,X(:,2));
% title('The other population over time');
% xlabel('time');
% ylabel('x_2(t)');
% nexttile
% plot(X(:,1),X(:,2));
% title('Prey - Predator population cycle')
% xlabel('x_1(t)') 
% ylabel('x_2(t)') 

%% ExplicitEuler
clear
clc

a = 1;
b = 1;
x0 = [2; 2];
N = 10^5; % Number of steps
options = odeset('RelTol',1.0e-6,'AbsTol',1.0e-6);
[T,X]=ExplicitEulerFixedStepSize(@PreyPredator,0, 50, N, x0, a,b );

figure;
figs3 = tiledlayout(2,2);
figs3.Padding = 'compact';
figs3.TileSpacing = 'compact';
title(figs3,'Prey-Predator solved with ExplicitEuler');
nexttile                % First subplot
plot(T,X(:,1));
title('One population over time');
xlabel('time') 
ylabel('x_1(t)') 
nexttile
plot(T,X(:,2));
title('The other population over time');
xlabel('time');
ylabel('x_2(t)');
nexttile
plot(X(:,1),X(:,2));
title('Prey - Predator population cycle')
xlabel('x_1(t)') 
ylabel('x_2(t)')

%% μ small ??, ImplicitEulerFixedStepSize.m + NewtonsMethodODE.m + VanderPolFunJac.m
clear
clc

a = 1;
b = 1;
x0 = [2; 2];
t0=0;
tN=50;
N = 10^5;            % Number of steps, with 10^3 steps it does not converge.   
options = odeset('RelTol',1.0e-6,'AbsTol',1.0e-6); % The options are not needed as input in the custom ExplicitEuler function
[T,X]=ImplicitEulerFixedStepSize(@PreyPredatorFunJac, [t0, tN], N, x0, a,b);

plotter_prey_pred(T,X,a,b,"ImplicitEuler")


