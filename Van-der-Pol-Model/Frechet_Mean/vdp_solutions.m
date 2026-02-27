function vdp_solutions()
% vdp_interactive  Interactive Van der Pol (slow-fast form) slider visualization
% Save this file as vdp_interactive.m and run by typing: vdp_interactive
clear; close all; clc;

set(groot, ...
    'DefaultTextInterpreter','latex', ...
    'DefaultAxesTickLabelInterpreter','latex', ...
    'DefaultLegendInterpreter','latex', ...
    'DefaultAxesFontName','Times', ...
    'DefaultTextFontName','Times', ...
    'DefaultAxesFontSize', 12, ...
    'DefaultTextFontSize', 12, ...
    'DefaultLegendFontSize', 11, ...
    'DefaultLineLineWidth', 1.2, ...
    'DefaultFigureColor','w');

% Parameters
epsilon = 0.1;                      % slow timescale
a = -0.5;
x0 = [-1; 0];                        % initial condition
tspan = [0 500];                     % long enough to reach steady-state

% Precompute trajectories for all a-values (for smooth slider)

vdp = @(t, X) [X(2) - X(1)^3/3 + X(1);
               epsilon * (a - X(1))];
[t, X] = ode15s(vdp, tspan, x0);

% ---- Remove transient ----
t_cut = 200;                       % discard first 200 time units
idx = t > t_cut;

t_ss = t(idx);
x_ss = X(idx,1);
y_ss = X(idx,2);


% ---- Create Figure ----
fig = figure('Name','Van der Pol Solutions','NumberTitle','off','Units','normalized','Position',[0.2 0.2 0.6 0.6]);
ax = axes('Parent', fig);
hold(ax, 'on');

axis(ax, [-2.5 2.5 -2 2]);
grid on;

vdp_curve = plot(ax, x_ss, y_ss, 'LineWidth', 2, 'DisplayName', 'Limit cycle');

% ---- Add four clean arrow heads to indicate direction ----

nArrows = 4;

% pick four roughly equally spaced locations
idx_arrows = [1,17,35, 55];

pt1 = ax.Position;
xlim = ax.XLim;
ylim = ax.YLim;

arrow_color = vdp_curve.Color;
for k = 1:nArrows
    
    i  = idx_arrows(k);
    i2 = i + 3;  % small forward offset
    
    % Convert to normalized figure coordinates
    x1 = (x_ss(i)  - xlim(1)) / diff(xlim) * pt1(3) + pt1(1);
    y1 = (y_ss(i)  - ylim(1)) / diff(ylim) * pt1(4) + pt1(2);

    x2 = (x_ss(i2) - xlim(1)) / diff(xlim) * pt1(3) + pt1(1);
    y2 = (y_ss(i2) - ylim(1)) / diff(ylim) * pt1(4) + pt1(2);

    annotation(fig,'arrow',[x1 x2],[y1 y2],'Color',arrow_color,'LineWidth',2);
end


% Plot static nullcline (x' = 0)
xgrid = linspace(-3,3,400);
y_nullcline = xgrid.^3/3 - xgrid;
plot(ax, xgrid, y_nullcline, 'k--', 'LineWidth', 1.5, 'DisplayName', 'x''=0');

exportgraphics(gcf, 'Figures_vdp/vdp_system.pdf', ...
        'ContentType','vector');

end