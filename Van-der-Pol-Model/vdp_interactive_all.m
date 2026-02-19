function vdp_interactive_all()
% vdp_interactive_all  Interactive visualization for a generalized slow–fast system:
%
%   x' = b*y - c*x^3/3 + d*x + e
%   y' = ε*(a - x) + f
%
% Sliders for a, b, c, d, e, f ∈ [-1, 3],  ε ∈ [0.01, 2].
% Run by typing:  vdp_interactive_all

clear; close all; clc;

% --- Parameter ranges ---
param_range = [-1, 3];
epsilon_range = [0.01, 2];

% --- Initial conditions ---
x0 = [-1; 0];
tspan = [0 300];

% --- Create figure and axis ---
fig = figure('Name','Generalized Van der Pol (Interactive, with f)', ...
    'NumberTitle','off','Units','normalized','Position',[0.05 0.05 0.9 0.9]);
ax = axes('Parent', fig);
hold(ax, 'on');
xlabel(ax, 'x'); ylabel(ax, 'y');
axis(ax, [-3 3 -4 4]);
grid on;
title(ax, 'x'' = b y - c x^3/3 + d x + e,   y'' = ε (a - x) + f');

% --- Adjust axes position to make room for sliders on the right ---
ax.Position = [0.05 0.05 0.65 0.9];  % left, bottom, width, height

% --- Vector field grid ---
[xv, yv] = meshgrid(linspace(-3,3,20), linspace(-4,4,20));

% --- Initial parameter values ---
a = 0; b = 1; c = 1; d = 1; e = 0; f = 0; epsilon = 0.1;

% --- Initial vector field ---
u = b*yv - c*xv.^3/3 + d*xv + e;
v = epsilon*(a - xv) + f;
m = sqrt(u.^2 + v.^2); u = u./m; v = v./m;
h_quiver = quiver(ax, xv, yv, u, v, 0.5, 'Color',[0.7 0.7 0.7],'HandleVisibility','off');

% --- Static x'-nullcline (depends on parameters dynamically) ---
xgrid = linspace(-3,3,400);
y_nullcline = (c*xgrid.^3/3 - d*xgrid - e)/b;  % from x'=0
h_xnull = plot(ax, xgrid, y_nullcline, 'k--', 'LineWidth', 1.5, 'DisplayName','x''=0');

% --- Dynamic elements ---
h_yline = xline(ax, a + f/epsilon, 'b--', 'LineWidth', 1.5, 'DisplayName','y''=0');
h_traj = plot(ax, NaN, NaN, 'r-', 'LineWidth', 2, 'DisplayName','Trajectory');
h_eq = plot(ax, NaN, NaN, 'ko', 'MarkerFaceColor','k', 'DisplayName','Equilibrium');
legend(ax, 'Location','best');

% --- Helper function to create sliders + labels ---
function [slider, label] = make_slider(ypos, range, init_val, name, slider_left, slider_width, label_left, label_width)
    % Create slider
    slider = uicontrol('Parent', fig, 'Style', 'slider', ...
        'Units', 'normalized', ...
        'Position', [slider_left ypos slider_width 0.04], ...
        'Min', range(1), 'Max', range(2), ...
        'Value', init_val, ...
        'SliderStep', [1/80 0.1]);
    
    % Create label above the slider
    label = uicontrol('Style', 'text', 'Units', 'normalized', ...
        'Position', [label_left ypos+0.03 label_width 0.04], ...  % slightly above
        'String', sprintf('%s = %.2f', name, init_val), ...
        'FontSize', 12, 'HorizontalAlignment', 'left');
end

% --- Slider layout on right side ---
slider_width = 0.2;
label_width  = 0.2;
slider_left  = 0.72;
label_left   = slider_left + slider_width + 0.02;

[y_eps, y_f, y_e, y_d, y_c, y_b, y_a] = deal(0.02, 0.08, 0.14, 0.20, 0.26, 0.32, 0.38);

[slider_a, label_a] = make_slider(y_a, param_range, a, 'a', slider_left, slider_width, label_left, label_width);
[slider_b, label_b] = make_slider(y_b, param_range, b, 'b', slider_left, slider_width, label_left, label_width);
[slider_c, label_c] = make_slider(y_c, param_range, c, 'c', slider_left, slider_width, label_left, label_width);
[slider_d, label_d] = make_slider(y_d, param_range, d, 'd', slider_left, slider_width, label_left, label_width);
[slider_e, label_e] = make_slider(y_e, param_range, e, 'e', slider_left, slider_width, label_left, label_width);
[slider_f, label_f] = make_slider(y_f, param_range, f, 'f', slider_left, slider_width, label_left, label_width);

slider_eps = uicontrol('Parent', fig, 'Style', 'slider', ...
    'Units', 'normalized', 'Position', [slider_left y_eps slider_width 0.04], ...
    'Min', epsilon_range(1), 'Max', epsilon_range(2), ...
    'Value', epsilon, 'SliderStep', [1/66 0.1]);
label_eps = uicontrol('Style', 'text', 'Units', 'normalized', ...
    'Position', [label_left y_eps label_width 0.04], ...
    'String', sprintf('ε = %.2f', epsilon), 'FontSize', 12, 'HorizontalAlignment','left');

% --- Update Function ---
function updatePlot(~, ~)
    % Read slider values
    a = slider_a.Value;
    b = slider_b.Value;
    c = slider_c.Value;
    d = slider_d.Value;
    e = slider_e.Value;
    f = slider_f.Value;
    epsilon = slider_eps.Value;

    % Update vector field
    u = b*yv - c*xv.^3/3 + d*xv + e;
    v = epsilon*(a - xv) + f;
    m = sqrt(u.^2 + v.^2); u = u./m; v = v./m;
    set(h_quiver, 'UData', u, 'VData', v);

    % Update x'-nullcline
    y_null = (c*xgrid.^3/3 - d*xgrid - e)/b;
    set(h_xnull, 'YData', y_null);

    % Update y'-nullcline
    h_yline.Value = a + f/epsilon;

    % Update equilibrium
    x_eq = a + f/epsilon;
    y_eq = (c*x_eq^3/3 - d*x_eq - e)/b;
    set(h_eq, 'XData', x_eq, 'YData', y_eq);

    % Integrate trajectory
    vdp = @(t,X) [b*X(2) - c*X(1)^3/3 + d*X(1) + e;
                  epsilon*(a - X(1)) + f];
    [t, X] = ode15s(vdp, tspan, x0);
    set(h_traj, 'XData', X(:,1), 'YData', X(:,2));

    % Update labels with variable names
    label_a.String = sprintf('a = %.2f', a);
    label_b.String = sprintf('b = %.2f', b);
    label_c.String = sprintf('c = %.2f', c);
    label_d.String = sprintf('d = %.2f', d);
    label_e.String = sprintf('e = %.2f', e);
    label_f.String = sprintf('f = %.2f', f);
    label_eps.String = sprintf('ε = %.2f', epsilon);

    drawnow;
end

% --- Assign Callbacks ---
sliders = [slider_a, slider_b, slider_c, slider_d, slider_e, slider_f, slider_eps];
for s = sliders
    s.Callback = @updatePlot;
end

% --- Initialize Plot ---
updatePlot();

% --- Reset Button ---
reset_button = uicontrol('Parent', fig, 'Style', 'pushbutton', ...
    'Units', 'normalized', ...
    'Position', [slider_left 0.45 0.2 0.05], ...
    'String', 'Reset', 'FontSize', 12, ...
    'Callback', @resetAll);

function resetAll(~,~)
    % Reset sliders to initial values
    slider_a.Value   = 0;
    slider_b.Value   = 1;
    slider_c.Value   = 1;
    slider_d.Value   = 1;
    slider_e.Value   = 0;
    slider_f.Value   = 0;
    slider_eps.Value = 0.1;

    % Force MATLAB to update the slider thumb positions
    drawnow;   % <--- THIS is the missing part

    % Force sliders to update visually by firing their callback
    updatePlot(slider_a, []);

end

end
