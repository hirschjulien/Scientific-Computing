 function [fig1,fig2,fig3] = plotter_cstr(T,X,p,solver)
% WORK IN PROGRESS

filename = 'vdp_mu'+ string(mu)+'_'+solver+'.pdf';

% Turn grid on for all axes
set(groot, 'defaultAxesXGrid', 'on');
set(groot, 'defaultAxesYGrid', 'on');

% Set default font size for axes labels and titles
set(groot, 'defaultAxesFontSize', 18);
set(groot, 'defaultAxesTitleFontSizeMultiplier', 1.3);
set(groot, 'defaultAxesLabelFontSizeMultiplier', 1.3);

% (Optional) Set default line width and font name
set(groot, 'defaultLineLineWidth', 1.6);
% set(groot, 'defaultAxesFontName', 'Arial');

set(groot, 'defaultLineLineWidth', 1.6);

figs1 = tiledlayout(3,1);
figs1.Padding = 'compact';
figs1.TileSpacing = 'compact';
ax1 = nexttile;
plot(T,X(:,1));
% pbaspect(ax1, [1 1 1])
% title('Solved with ode45');
ylabel('$$C_A(t)$$','interpreter','latex')
ax2 = nexttile;
plot(T,X(:,2));
% pbaspect(ax2,[1 1 1])
% title('Solved with ode45');
% xlabel('time');
ylabel('$$C_B(t)$$','interpreter','latex');
ax3 = nexttile;
plot(T,X(:,3));
ylabel('$$T(t)$$','interpreter','latex');
xlabel('$$t [min]$$', 'interpreter', 'latex') 


% COPYPASTA BELOW

fig1 = figure(1);
plot(T,X(:,1));
pbaspect([1 1 1])
xlim([min(T), max(T)])
% title('Solved with ode45');
xlabel('time') 
ylabel('$$x_1(t)$$','interpreter','latex')
exportgraphics(fig1,'x1_t_'+filename,'ContentType','vector')


fig2 = figure(2);
plot(T,X(:,2));
xlim([min(T), max(T)])
pbaspect([1 1 1])
% pbaspect(fig2,[1 1 1])
% title('Solved with ode45');
xlabel('time');
ylabel('$$x_2(t)$$','interpreter','latex');
exportgraphics(fig2,'x2_t_'+filename,'ContentType','vector')

fig3 = figure(3);
plot(X(:,1),X(:,2));
pbaspect([1 1 1])
% pbaspect(fig3, [1 1 1])
% title('VDP - Solved with ode45')
xlabel('$$x_1(t)$$','interpreter','latex')
ylabel('$$x_2(t)$$','interpreter','latex')
exportgraphics(fig3,'x2_x1_'+filename','ContentType','vector')

disp("VdP Plotter finished!")
