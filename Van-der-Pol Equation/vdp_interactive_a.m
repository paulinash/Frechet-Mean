function vdp_interactive_a()
% vdp_interactive  Interactive Van der Pol (slow-fast form) slider visualization
% Save this file as vdp_interactive.m and run by typing: vdp_interactive
clear; close all; clc;

% Parameters
epsilon = 0.1;                      % slow timescale
a_values = linspace(-1.5, 1.5, 91);  % values of 'a' to explore
x0 = [-1; 0];                        % initial condition
tspan = [0 500];                     % long enough to reach steady-state

% Precompute trajectories for all a-values (for smooth slider)
trajectories = cell(length(a_values), 1);
for i = 1:length(a_values)
    a = a_values(i);
    vdp = @(t, X) [X(2) - X(1)^3/3 + X(1);
                   epsilon * (a - X(1))];
    [t, X] = ode15s(vdp, tspan, x0);
    trajectories{i}.t = t;
    trajectories{i}.x = X(:,1);
    trajectories{i}.y = X(:,2);
end

% ---- Create Figure ----
fig = figure('Name','Van der Pol Interactive','NumberTitle','off','Units','normalized','Position',[0.2 0.2 0.6 0.6]);
ax = axes('Parent', fig);
hold(ax, 'on');
xlabel(ax, 'x');
ylabel(ax, 'y');
title(ax, 'Phase Portrait with Nullclines and Trajectory');
axis(ax, [-3 3 -4 4]);
grid on;


% Grid for vector field
[xv, yv] = meshgrid(linspace(-3, 3, 20), linspace(-4, 4, 20));

% --- Initial vector field (for a = a_values(1)) ---
a0 = a_values(1);
u = yv - xv.^3/3 + xv;
v = epsilon * (a0 - xv);
m = sqrt(u.^2 + v.^2);
u = u ./ m; v = v ./ m; % normalize
h_quiver = quiver(ax, xv, yv, u, v, 0.5, 'Color', [0.7 0.7 0.7], ...
                  'HandleVisibility','off');


% Plot static nullcline (x' = 0)
xgrid = linspace(-3,3,400);
y_nullcline = xgrid.^3/3 - xgrid;
plot(ax, xgrid, y_nullcline, 'k--', 'LineWidth', 1.5, 'DisplayName', 'x''=0');



% Add dynamic elements (that change with slider)
h_yline = xline(ax, a_values(1), 'b--', 'LineWidth', 1.5, 'DisplayName', 'y''=0');
h_traj   = plot(ax, NaN, NaN, 'r-', 'LineWidth', 2, 'DisplayName', 'Trajectory');
h_eq     = plot(ax, NaN, NaN, 'ko', 'MarkerFaceColor','k', 'DisplayName','Equilibrium');

legend(ax, 'Location', 'best');

% ---- Create Slider ----
slider = uicontrol('Parent', fig, 'Style', 'slider', ...
    'Units', 'normalized', ...
    'Position', [0.2 0.02 0.6 0.05], ...
    'Min', 1, 'Max', length(a_values), ...
    'Value', 1, 'SliderStep', [1/(length(a_values)-1) 0.1]);

% Text label for current 'a'
label = uicontrol('Style', 'text', 'Units', 'normalized', ...
    'Position', [0.82 0.02 0.15 0.05], ...
    'String', sprintf('a = %.2f', a_values(1)), ...
    'FontSize', 12);

% ---- Nested callback (has access to a_values, trajectories, handles) ----
    function updatePlot(val)
        idx = round(val);
        % clamp index
        idx = max(1, min(length(a_values), idx));
        a = a_values(idx);

        %%% >>> NEW <<< : Update the vector field dynamically
        u = yv - xv.^3/3 + xv;        % x' (same)
        v = epsilon * (a - xv);       % y' (changes with a)
        m = sqrt(u.^2 + v.^2);
        u = u ./ m; v = v ./ m;
        set(h_quiver, 'UData', u, 'VData', v);


     
        % Update dynamic elements
        h_yline.Value = a;                        % move vertical nullcline
          x_eq = a;
        y_eq = x_eq^3/3 - x_eq;
        set(h_traj, 'XData', trajectories{idx}.x, ...
                    'YData', trajectories{idx}.y);
        set(h_eq,   'XData', x_eq, 'YData', y_eq);
        label.String = sprintf('a = %.2f', a);

        drawnow;
    end

% Assign callback to slider (use Value passed to nested updatePlot)
slider.Callback = @(src, ~) updatePlot(src.Value);

% Initialize display
updatePlot(1);

end