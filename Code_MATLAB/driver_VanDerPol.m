%% Non-stiff problem (µ small) - explicit solver - ode45
clear
clc
mu = 10;
x0 = [2.0; 0.0];
options = odeset('RelTol',1.0e-6,'AbsTol',1.0e-6);
[T,X]=ode45(@VanDerPol,[0 5*mu],x0,options,mu);

p.mu=mu;
p.solver="ode45";
p.export_gfx=false;
p.abstol(1)=1.0e-6;
plotter_vdp(T,X,p);

% figure;
% figs1 = tiledlayout(1,3);
% figs1.Padding = 'compact';
% figs1.TileSpacing = 'compact';
% ax1 = nexttile;
% plot(T,X(:,1));
% pbaspect(ax1, [1 1 1])
% title('Solved with ode45');
% xlabel('time') 
% ylabel('$$x_1(t)$$','interpreter','latex')
% ax2 = nexttile;
% plot(T,X(:,2));
% pbaspect(ax2,[1 1 1])
% title('Solved with ode45');
% xlabel('time');
% ylabel('$$x_2(t)$$','interpreter','latex');
% ax3 = nexttile;
% plot(X(:,1),X(:,2));
% pbaspect(ax3, [1 1 1])
% title('VDP -Solved with ode45')
% xlabel('$$x_1(t)$$','interpreter','latex')
% ylabel('$$x_2(t)$$','interpreter','latex')
% % set(figs1,'PaperPositionMode','auto');
% %saveas(figs1,'vdp_mu3_ode45.pdf')
% exportgraphics(figs1,'vdp_mu3_ode45.pdf','ContentType','vector')

%% Stiff problem (µ ~large) - implicit solver - ode15s
clear
clc
mu = 10;          % Ode four five breaks at μ=1000 and takes too long
x0 = [2.0; 0.0];
options = odeset('Jacobian',@VanDerPolJac,'RelTol',1.0e-6,'AbsTol',1.0e-6);
[T,X]=ode15s(@VanDerPol,[0 5*mu],x0,options,mu);

p.mu=mu;
p.solver="ode15s";
p.export_gfx=false;
p.abstol_sing=1.0e-6;
plotter_vdp(T,X,p);

%% Stiff problem μ10, ode15s - Triple comparison
clear
clc
mu = 10;          % Ode four five breaks at μ=1000 and takes too long
x0 = [2.0; 0.0];
p.abstol(1)=1.0e-3;
p.abstol(2)=1.0e-6;
p.abstol(3)=1.0e-9;
options1 = odeset('Jacobian',@VanDerPolJac,'RelTol',p.abstol(1),'AbsTol',p.abstol(1));
[T1,X1]=ode15s(@VanDerPol,[0 5*mu],x0,options1,mu);
options2 = odeset('Jacobian',@VanDerPolJac,'RelTol',p.abstol(2),'AbsTol',p.abstol(2));
[T2,X2]=ode15s(@VanDerPol,[0 5*mu],x0,options2,mu);
options3 = odeset('Jacobian',@VanDerPolJac,'RelTol',p.abstol(3),'AbsTol',p.abstol(3));
[T3,X3]=ode15s(@VanDerPol,[0 5*mu],x0,options3,mu);

p.mu=mu;
p.solver="ode15s";
p.export_gfx=true;

ttr=cell(3,1);
xtr=cell(3,1);
ttr{1}=T1;
ttr{2}=T2;
ttr{3}=T3;
xtr{1}=X1;
xtr{2}=X2;
xtr{3}=X3;
triplotter_vdp(ttr,xtr,p);


%% μ small, ExplicitEulerFixedStepSize.m
clear
clc

mu = 10;        % In some slides μ small is μ=3!
x0 = [2.0; 0.0];
N = 10^4;            % Number of steps, with 10^3 steps it does not converge.   

[T,X]=ExplicitEulerFixedStepSize(@VanDerPol, [0, 5*mu], N, x0, mu);

% PLOTS
figure;
figs3 = tiledlayout(2,2);
figs3.Padding = 'compact';
figs3.TileSpacing = 'compact';
title(figs3, sprintf('Solver: ExplicitEuler, μ=%.1f', mu'));
nexttile                % First subfigure
plot(T,X(:,1));
%title('Solved with ExplicitEuler');
xlabel('time') 
ylabel('x_1(t)')
xlim([0,5*mu]);
nexttile                % Second subfigure
plot(T,X(:,2));
%title('Solved with ExplicitEuler');
xlabel('time');
ylabel('x_2(t)');
xlim([0,5*mu]);
nexttile                % Third subfigure
plot(X(:,1),X(:,2));
%title('VDP -Solved with ExplicitEuler')
xlabel('x_1(t)') 
ylabel('x_2(t)')

%% μ small ??, ImplicitEulerFixedStepSize.m + NewtonsMethodODE.m + VanderPolFunJac.m
clear
clc

mu = 10;
t0 = 0;
tf = 5*mu;
x0 = [2.0; 0.0];
N = 10^4;            % Number of steps, with 10^3 steps it does not converge.   
%options = odeset('RelTol',1.0e-6,'AbsTol',1.0e-6); % The options are not needed as input in the custom ExplicitEuler function
[T,X]=ImplicitEulerFixedStepSize(@VanDerPolFunJac, [t0, tf], N, 1.e-8, x0, mu);


% PLOTS
figure;
figs4 = tiledlayout(2,2);
figs4.Padding = 'compact';
figs4.TileSpacing = 'compact';
title(figs4, sprintf('Solver: ImplicitEuler, μ=%.1f', mu'));
nexttile                % First subfigure
plot(T,X(:,1));
%title('Solved with ImplicitEuler');
xlabel('time');
ylabel('x_1(t)');
xlim([0,5*mu]);
nexttile                % Second subfigure
plot(T,X(:,2));
%title('Solved with ExplicitEuler');
xlabel('time');
ylabel('x_2(t)');
xlim([0,5*mu]);
nexttile                % Third subfigure
plot(X(:,1),X(:,2));
%title('VDP -Solved with ExplicitEuler')
xlabel('x_1(t)') ;
ylabel('x_2(t)');
%% μ=10, ImplicitEulerFixedStepSize.m + NewtonsMethodODE.m + VanderPolFunJac.m - 3 Tolerances
clear
clc
mu = 10;          % Ode four five breaks at μ=1000 and takes too long
x0 = [2.0; 0.0];
N = 10^4; 
p.abstol(1)=1.0e-3;
p.abstol(2)=1.0e-4;
p.abstol(3)=1.0e-6;

[T1,X1]=ImplicitEulerFixedStepSize(@VanDerPolFunJac,[0 5*mu], N, p.abstol(1), x0, mu );

[T2,X2]=ImplicitEulerFixedStepSize(@VanDerPolFunJac,[0 5*mu], N, p.abstol(2), x0, mu );

[T3,X3]=ImplicitEulerFixedStepSize(@VanDerPolFunJac,[0 5*mu], N, p.abstol(3), x0, mu );

p.mu=mu;
p.solver="Impl.Euler_mid";
p.export_gfx=true;

ttr=cell(3,1);
xtr=cell(3,1);
ttr{1}=T1;
ttr{2}=T2;
ttr{3}=T3;
xtr{1}=X1;
xtr{2}=X2;
xtr{3}=X3;
triplotter_vdp_tols(ttr,xtr,p);

%% μ=10, ImplicitEulerFixedStepSize.m + NewtonsMethodODE.m + VanderPolFunJac.m - 3 Steps
clear
clc
mu = 10;          % Ode four five breaks at μ=1000 and takes too long
x0 = [2.0; 0.0];
p.N(1) = 1e4; 
p.N(2) = 1e5;
p.N(3) = 1e6;
tol = 1e6;

[T1,X1]=ImplicitEulerFixedStepSize(@VanDerPolFunJac,[0 5*mu], p.N(1), tol, x0, mu );

[T2,X2]=ImplicitEulerFixedStepSize(@VanDerPolFunJac,[0 5*mu], p.N(2), tol, x0, mu );

[T3,X3]=ImplicitEulerFixedStepSize(@VanDerPolFunJac,[0 5*mu], p.N(3), tol, x0, mu );

p.mu=mu;
p.solver="ImplEuler";
p.export_gfx=true;

ttr=cell(3,1);
xtr=cell(3,1);
ttr{1}=T1;
ttr{2}=T2;
ttr{3}=T3;
xtr{1}=X1;
xtr{2}=X2;
xtr{3}=X3;
triplotter_vdp_steps(ttr,xtr,p);
%% ADAPTIVE, μ=10, ExplicitEulerAdaptiveStep.m
% The adaptive step requires the following ADDITIONAL inputs:
% h0 = initial step guess
% abstol = Absolute tolerance
% reltol = Relative tolerance

clear
clc

mu = 3;
t0 = 0;
tf = 5*mu;
x0 = [2.0; 0.0];
% VAriables needed for adaptive step method:
abstol = 1e-6; % Usual / default value 
reltol = 1e-6;
% INITIAL STEP SIZE CALCULATION
h0 = InitialStepSize(@VanDerPol,t0,x0,mu);
% d0 = norm(x0);
% d1 = norm(feval(@VanDerPol,t0,x0,mu));
% if 
%     h0 = = .01(d0/d1);
[T,X]=ExplicitEulerAdaptiveStep(@VanDerPol, [t0, tf], x0, h0, abstol, reltol, mu);

% PLOTS
figure;
figs4 = tiledlayout(2,2);
figs4.Padding = 'compact';
figs4.TileSpacing = 'compact';
title(figs4, sprintf('Solver: AdaptiveExplicitEuler, μ=%.1f', mu'));
nexttile                % First subfigure
plot(T,X(:,1));
xlabel('time');
ylabel('x_1(t)');
xlim([0,5*mu]);
nexttile                % Second subfigure
plot(T,X(:,2));
xlabel('time');
ylabel('x_2(t)');
xlim([0,5*mu]);
nexttile                % Third subfigure
plot(X(:,1),X(:,2));
xlabel('x_1(t)') ;
ylabel('x_2(t)');

%% RK4 Fixed, μ=10, CLassicalRungeKuttaFixedStep.m
clear
clc

mu = 10;
t0 = 0;
tf = 5*mu;
x0 = [2.0; 0.0];
N=1e4;

[T,X]=ClassicalRungeKuttaFixedStep(@VanDerPol, [t0, tf], x0, N, mu);

% PLOTS
figure;
figs4 = tiledlayout(2,2);
figs4.Padding = 'compact';
figs4.TileSpacing = 'compact';
title(figs4, sprintf('Solver: Adaptive Classical RK, μ=%.1f', mu'));
nexttile                % First subfigure
plot(T,X(:,1));
xlabel('time');
ylabel('x_1(t)');
xlim([0,5*mu]);
nexttile                % Second subfigure
plot(T,X(:,2));
xlabel('time');
ylabel('x_2(t)');
xlim([0,5*mu]);
nexttile                % Third subfigure
plot(X(:,1),X(:,2));
xlabel('x_1(t)') ;
ylabel('x_2(t)');


%% ADAPTIVE, μ=10, CLassicalRungeKuttaAdaptiveStep.m
% The adaptive step requires the following ADDITIONAL inputs:
% h0 = initial step guess
% abstol = Absolute tolerance
% reltol = Relative tolerance

clear
clc

mu = 10;
t0 = 0;
tf = 5*mu;
x0 = [2.0; 0.0];
% VAriables needed for adaptive step method:
abstol = 1e-6; % Usual / default value is 1e-6
reltol = 1e-6; % Which is 
% INITIAL STEP SIZE CALCULATION
h0 = InitialStepSize(@VanDerPol,t0,x0,mu);
% d0 = norm(x0);
% d1 = norm(feval(@VanDerPol,t0,x0,mu));
% if 
%     h0 = = .01(d0/d1);

[T,X]=ClassicalRungeKuttaAdaptiveStep(@VanDerPol, [t0, tf], x0, h0, abstol, reltol, mu);

% PLOTS
figure;
figs4 = tiledlayout(2,2);
figs4.Padding = 'compact';
figs4.TileSpacing = 'compact';
title(figs4, sprintf('Solver: Adaptive Classical RK, μ=%.1f', mu'));
nexttile                % First subfigure
plot(T,X(:,1));
xlabel('time');
ylabel('x_1(t)');
xlim([0,5*mu]);
nexttile                % Second subfigure
plot(T,X(:,2));
xlabel('time');
ylabel('x_2(t)');
xlim([0,5*mu]);
nexttile                % Third subfigure
plot(X(:,1),X(:,2));
xlabel('x_1(t)') ;
ylabel('x_2(t)');

%% ADAPTIVE, μ=10, ImplicitEulerAdaptiveStep.m
% The adaptive step requires the following ADDITIONAL inputs:
% h0 = initial step guess
% abstol = Absolute tolerance
% reltol = Relative tolerance

clear
clc

mu = 10;
t0 = 0;
tf = 5*mu;
x0 = [2.0; 0.0];
% VAriables needed for adaptive step method:
abstol = 1e-6; % Usual / default value is 1e-6
reltol = 1e-6; % Which is 
% INITIAL STEP SIZE CALCULATION
h0 = InitialStepSize(@VanDerPol,t0,x0,mu);
% d0 = norm(x0);
% d1 = norm(feval(@VanDerPol,t0,x0,mu));
% if 
%     h0 = = .01(d0/d1);

options = odeset('RelTol',1.0e-6,'AbsTol',1.0e-6); % The options are not needed as input in the custom ExplicitEuler function
[T,X]=ImplicitEulerAdaptiveStep(@VanDerPolFunJac, [t0, tf], x0, h0, abstol, reltol, mu);

% PLOTS
figure;
figs4 = tiledlayout(2,2);
figs4.Padding = 'compact';
figs4.TileSpacing = 'compact';
title(figs4, sprintf('Solver: Adaptive ImplEul, μ=%.1f', mu'));
nexttile                % First subfigure
plot(T,X(:,1));
xlabel('time');
ylabel('x_1(t)');
xlim([0,5*mu]);
nexttile                % Second subfigure
plot(T,X(:,2));
xlabel('time');
ylabel('x_2(t)');
xlim([0,5*mu]);
nexttile                % Third subfigure
plot(X(:,1),X(:,2));
xlabel('x_1(t)') ;
ylabel('x_2(t)');


%% ADAPTIVE, COMPARISON μ=10, ImplicitEulerAdaptiveStep.m + NewtonsMethodODE.m + VanderPolFunJac.m - 3 Absolute Tolerances
clear
clc
mu = 10;
x0 = [2.0; 0.0];
t0 = 0;
tf = 5*mu;
% INITIAL STEP SIZE CALCULATION
h0 = InitialStepSize(@VanDerPol,t0,x0,mu);
% Solver parameters
reltol=1e-3;
p.abstol(1)=1.0e-3;
p.abstol(2)=1.0e-4;
p.abstol(3)=1.0e-6;

%ImplicitEulerAdaptiveStep(FunJac, tspan, x0, h0, abstol, reltol, varargin)
[T1,X1]=ImplicitEulerAdaptiveStep(@VanDerPolFunJac,[t0 tf], x0, h0, p.abstol(1), reltol, mu);

[T2,X2]=ImplicitEulerAdaptiveStep(@VanDerPolFunJac,[t0 tf], x0, h0, p.abstol(2), reltol, mu);

[T3,X3]=ImplicitEulerAdaptiveStep(@VanDerPolFunJac,[t0 tf], x0, h0, p.abstol(3), reltol, mu);

p.mu=mu;
p.solver="Impl.EulerAdaptive_mid";
p.export_gfx=true;

ttr=cell(3,1);
xtr=cell(3,1);
ttr{1}=T1;
ttr{2}=T2;
ttr{3}=T3;
xtr{1}=X1;
xtr{2}=X2;
xtr{3}=X3;
triplotter_vdp_tols(ttr,xtr,p);

%% ADAPTIVE, COMPARISON μ=10, ImplicitEulerAdaptiveStep.m + NewtonsMethodODE.m + VanderPolFunJac.m - 3 RELATIVE Tolerances
clear
clc
mu = 10;
x0 = [2.0; 0.0];
t0 = 0;
tf = 5*mu;
% INITIAL STEP SIZE CALCULATION
h0 = InitialStepSize(@VanDerPol,t0,x0,mu);
% Solver parameters
abstol=1e-6;
p.reltol(1)=1.0e-3;
p.reltol(2)=1.0e-4;
p.reltol(3)=1.0e-6;

%ImplicitEulerAdaptiveStep(FunJac, tspan, x0, h0, abstol, reltol, varargin)
[T1,X1]=ImplicitEulerAdaptiveStep(@VanDerPolFunJac,[t0 tf], x0, h0, abstol, p.reltol(1), mu);

[T2,X2]=ImplicitEulerAdaptiveStep(@VanDerPolFunJac,[t0 tf], x0, h0, abstol, p.reltol(2), mu);

[T3,X3]=ImplicitEulerAdaptiveStep(@VanDerPolFunJac,[t0 tf], x0, h0, abstol, p.reltol(3), mu);

p.mu=mu;
p.solver="Impl.EulerAdaptive_MidRelTols";
p.export_gfx=true;

ttr=cell(3,1);
xtr=cell(3,1);
ttr{1}=T1;
ttr{2}=T2;
ttr{3}=T3;
xtr{1}=X1;
xtr{2}=X2;
xtr{3}=X3;
triplotter_vdp_reltols(ttr,xtr,p);

%% RK4 COMPARISONS
%% μ=10, RK4 Requires: ClassicalRungeKuttaFixedStep.m + ClasicalRungeKuttaStep.m - 3 Steps
clear
clc
mu = 10;          % Ode four five breaks at μ=1000 and takes too long
x0 = [2.0; 0.0];
p.N(1) = 1e4; 
p.N(2) = 1e5;
p.N(3) = 1e6;
tol = 1e6;

[T1,X1]=ClassicalRungeKuttaFixedStep(@VanDerPol,[0 5*mu], x0, p.N(1), mu );

[T2,X2]=ClassicalRungeKuttaFixedStep(@VanDerPol,[0 5*mu], x0, p.N(2), mu );

[T3,X3]=ClassicalRungeKuttaFixedStep(@VanDerPol,[0 5*mu], x0, p.N(3), mu );

p.mu=mu;
p.solver="RK4_456";
p.export_gfx=true;

ttr=cell(3,1);
xtr=cell(3,1);
ttr{1}=T1;
ttr{2}=T2;
ttr{3}=T3;
xtr{1}=X1;
xtr{2}=X2;
xtr{3}=X3;
triplotter_vdp_steps(ttr,xtr,p);

%% ADAPTIVE RK4 μ=10, ClassicalRungeKuttaAdaptiveStep.m + ClasicalRungeKuttaStep.m  - 3 Absolute Tolerances
clear
clc
mu = 10;
x0 = [2.0; 0.0];
t0 = 0;
tf = 5*mu;
% INITIAL STEP SIZE CALCULATION
h0 = InitialStepSize(@VanDerPol,t0,x0,mu);
% Solver parameters
reltol=1e-3;
p.abstol(1)=1.0e-9;
p.abstol(2)=1.0e-12;
p.abstol(3)=1.0e-14;

%ImplicitEulerAdaptiveStep(FunJac, tspan, x0, h0, abstol, reltol, varargin)
[T1,X1]=ClassicalRungeKuttaAdaptiveStep(@VanDerPol,[t0 tf], x0, h0, p.abstol(1), reltol, mu);

[T2,X2]=ClassicalRungeKuttaAdaptiveStep(@VanDerPol,[t0 tf], x0, h0, p.abstol(2), reltol, mu);

[T3,X3]=ClassicalRungeKuttaAdaptiveStep(@VanDerPol,[t0 tf], x0, h0, p.abstol(3), reltol, mu);

p.mu=mu;
p.solver="RK4_high";
p.export_gfx=true;

ttr=cell(3,1);
xtr=cell(3,1);
ttr{1}=T1;
ttr{2}=T2;
ttr{3}=T3;
xtr{1}=X1;
xtr{2}=X2;
xtr{3}=X3;
triplotter_vdp_tols(ttr,xtr,p);

%% ADAPTIVE RK4 μ=10, ClassicalRungeKuttaAdaptiveStep.m + ClasicalRungeKuttaStep.m - 3 RELATIVE Tolerances
clear
clc
mu = 10;
x0 = [2.0; 0.0];
t0 = 0;
tf = 5*mu;
% INITIAL STEP SIZE CALCULATION
h0 = InitialStepSize(@VanDerPol,t0,x0,mu);
% Solver parameters
abstol=1e-6;
p.reltol(1)=1.0e-8;
p.reltol(2)=1.0e-10;
p.reltol(3)=1.0e-12;

%ImplicitEulerAdaptiveStep(FunJac, tspan, x0, h0, abstol, reltol, varargin)
[T1,X1]=ImplicitEulerAdaptiveStep(@VanDerPolFunJac,[t0 tf], x0, h0, abstol, p.reltol(1), mu);

[T2,X2]=ImplicitEulerAdaptiveStep(@VanDerPolFunJac,[t0 tf], x0, h0, abstol, p.reltol(2), mu);

[T3,X3]=ImplicitEulerAdaptiveStep(@VanDerPolFunJac,[t0 tf], x0, h0, abstol, p.reltol(3), mu);

p.mu=mu;
p.solver="RK4_HighRelTols";
p.export_gfx=true;

ttr=cell(3,1);
xtr=cell(3,1);
ttr{1}=T1;
ttr{2}=T2;
ttr{3}=T3;
xtr{1}=X1;
xtr{2}=X2;
xtr{3}=X3;
triplotter_vdp_reltols(ttr,xtr,p);