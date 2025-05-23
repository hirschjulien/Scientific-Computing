function [fig1,fig2,fig3] = triplotter_vdp_steps(T,X,p)
T1 = T{1};
T2 = T{2};
T3 = T{3};
X1 = X{1};
X2 = X{2};
X3 = X{3};

mu = p.mu;
solver = p.solver;
export_gfx =p.export_gfx;
% tols = p.N;
filename = 'vdp_mu'+ string(mu)+'_'+solver+'_maxsteps.pdf';

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


fig1 = figure(1);
plot(nan, nan, 'w', "DisplayName", "N:")
hold on
plot(T1,X1(:,1), "DisplayName", sprintf('%.0e',p.N(1)) );
plot(T2,X2(:,1), "DisplayName", sprintf('%.0e',p.N(2)) );
plot(T3,X3(:,1), "DisplayName", sprintf('%.0e',p.N(3)) );
hold off
pbaspect([1 1 1])
legend( 'Interpreter', 'latex', 'Location', 'southwest')
xlim([min(T1), max(T2)])
xlabel('$$t$$','interpreter','latex');
ylabel('$$x_1(t)$$','interpreter','latex');
if export_gfx
    exportgraphics(fig1,'x1_t_'+filename,'ContentType','vector')
end


fig2 = figure(2);
plot(nan, nan, 'w', "DisplayName", "N:")
hold on
plot(T1,X1(:,2), "DisplayName", sprintf('%.0e',p.N(1)) );
plot(T2,X2(:,2), "DisplayName", sprintf('%.0e',p.N(2)) );
plot(T3,X3(:,2), "DisplayName", sprintf('%.0e',p.N(3)) );
hold off
pbaspect([1 1 1])
legend( 'Interpreter', 'latex','Location', 'southwest')
xlim([min(T1), max(T2)])
pbaspect([1 1 1])
xlabel('$$t$$','interpreter','latex');
ylabel('$$x_2(t)$$','interpreter','latex');
if export_gfx
    exportgraphics(fig2,'x2_t_'+filename,'ContentType','vector')
end

fig3 = figure(3);
plot(nan, nan, 'w', "DisplayName", "N:")
hold on
plot(X1(:,1),X1(:,2), "DisplayName", sprintf('%.0e',p.N(1)) );
plot(X2(:,1),X2(:,2), "DisplayName", sprintf('%.0e',p.N(2)) );
plot(X3(:,1),X3(:,2), "DisplayName", sprintf('%.0e',p.N(3)) );
hold off
pbaspect([1 1 1])
legend( 'Interpreter', 'latex','Location', 'northwest')
xlabel('$$x_1(t)$$','interpreter','latex')
ylabel('$$x_2(t)$$','interpreter','latex')
if export_gfx
    exportgraphics(fig3,'x2_x1_'+filename','ContentType','vector')
end

disp("VdP Plotter finished!")