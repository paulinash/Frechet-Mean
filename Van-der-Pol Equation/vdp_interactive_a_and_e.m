function vdp_interactive_a_and_e ()
% vdp_interactive  Interactive Van der Pol (slow-fast form) slider visualization
% Extended version: vary both a and epsilon interactively.

clear; close all; clc;

% Parameters
epsilon_values = linspace(0.01, 2.0, 67);   % includes 1.0 and 2.0 exactly
a_values = linspace(-1.5, 1.5, 91);           % range for a
x0 = [-1; 0];                                 % initial condition
tspan = [0 500];                              % long enough to reach steady state

% Precompute trajectories for all (a, epsilon) combinations
% To save time, we'll precompute for all 'a' and a *subset* of epsilon
% but you can increase density if you want.
n_eps = length(epsilon_values);
n_a   = length(a_values);
trajectories = cell(n_eps, n_a);

fprintf('Precomputing trajectories... this may take a minute.\n');
for ie = 1:n_eps
    epsilon = epsilon_values(ie);
    for ia = 1:n_a
        a = a_values(ia);
        vdp = @(t, X) [X(2) - X(1)^3/3 + X(1);
                       epsilon * (a - X(1))];
        [t, X] = ode15s(vdp, tspan, x0);
        trajectories{ie, ia}.t = t;
        trajectories{ie, ia}.x = X(:,1);
        trajectories{ie, ia}.y = X(:,2);
    end
end
fprintf('Done.\n');

% ---- Create Figure ----
fig = figure('Name','Van der Pol Interactive','NumberTitle','off','Units','normalized','Position',[0.15 0.15 0.7 0.7]);
ax = axes('Parent', fig);
hold(ax, 'on');
xlabel(ax, 'x');
ylabel(ax, 'y');
title(ax, 'Phase Portrait with Nullclines and Trajectory');
axis(ax, [-3 3 -4 4]);
grid on;

% Grid for vector field
[xv, yv] = meshgrid(linspace(-3, 3, 20), linspace(-4, 4, 20));

% --- Initial values ---
a_idx = 1;
eps_idx = 1;
a = a_values(a_idx);
epsilon = epsilon_values(eps_idx);

% Initial vector field
u = yv - xv.^3/3 + xv;
v = epsilon * (a - xv);
m = sqrt(u.^2 + v.^2);
u = u ./ m; v = v ./ m;
h_quiver = quiver(ax, xv, yv, u, v, 0.5, 'Color', [0.7 0.7 0.7], 'HandleVisibility','off');

% Plot static nullcline (x' = 0)
xgrid = linspace(-3,3,400);
y_nullcline = xgrid.^3/3 - xgrid;
plot(ax, xgrid, y_nullcline, 'k--', 'LineWidth', 1.5, 'DisplayName', 'x''=0');

% Dynamic elements
h_yline = xline(ax, a, 'b--', 'LineWidth', 1.5, 'DisplayName', 'y''=0');
h_traj   = plot(ax, NaN, NaN, 'r-', 'LineWidth', 2, 'DisplayName', 'Trajectory');
h_eq     = plot(ax, NaN, NaN, 'ko', 'MarkerFaceColor','k', 'DisplayName','Equilibrium');
legend(ax, 'Location', 'best');

% ---- Slider for 'a' ----
slider_a = uicontrol('Parent', fig, 'Style', 'slider', ...
    'Units', 'normalized', ...
    'Position', [0.2 0.02 0.6 0.04], ...
    'Min', 1, 'Max', n_a, ...
    'Value', a_idx, ...
    'SliderStep', [1/(n_a-1) 0.1]);

label_a = uicontrol('Style', 'text', 'Units', 'normalized', ...
    'Position', [0.82 0.02 0.15 0.04], ...
    'String', sprintf('a = %.2f', a_values(a_idx)), ...
    'FontSize', 12);

% ---- Slider for epsilon ----
slider_eps = uicontrol('Parent', fig, 'Style', 'slider', ...
    'Units', 'normalized', ...
    'Position', [0.2 0.08 0.6 0.04], ...
    'Min', 1, 'Max', n_eps, ...
    'Value', eps_idx, ...
    'SliderStep', [1/(n_eps-1) 0.1]);

label_eps = uicontrol('Style', 'text', 'Units', 'normalized', ...
    'Position', [0.82 0.08 0.15 0.04], ...
    'String', sprintf('ε = %.2f', epsilon_values(eps_idx)), ...
    'FontSize', 12);

% ---- Update Function ----
    function updatePlot(~, ~)
        a_idx = round(slider_a.Value);
        eps_idx = round(slider_eps.Value);
        a = a_values(a_idx);
        epsilon = epsilon_values(eps_idx);

        % Update vector field
        u = yv - xv.^3/3 + xv;
        v = epsilon * (a - xv);
        m = sqrt(u.^2 + v.^2);
        u = u ./ m; v = v ./ m;
        set(h_quiver, 'UData', u, 'VData', v);

        % Update trajectory and equilibrium
        x_eq = a;
        y_eq = x_eq^3/3 - x_eq;
        set(h_yline, 'Value', a);
        set(h_traj, 'XData', trajectories{eps_idx, a_idx}.x, ...
                    'YData', trajectories{eps_idx, a_idx}.y);
        set(h_eq, 'XData', x_eq, 'YData', y_eq);

        % Update labels
        label_a.String = sprintf('a = %.2f', a);
        label_eps.String = sprintf('ε = %.2f', epsilon);

        drawnow;
    end

% ---- Assign Callbacks ----
slider_a.Callback = @updatePlot;
slider_eps.Callback = @updatePlot;

% Initialize
updatePlot();

end