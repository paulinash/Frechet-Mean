% Van der Pol–type system (slow-fast form)
clear; close all; clc;

% Parameters
epsilon = 0.05;           % small parameter
a_values = [-1.5,-1,0,0.5,1,1.5]; % values of 'a' to test

% Create a colormap with the same number of colors as a_list
cmap = parula(length(a_values));  % you can also use jet, turbo, etc.

% Integration setup
tspan = [0 500];          % integrate long enough to reach steady behavior
x0 = [-2; 1];              % initial conditions


% Grid for nullclines
xgrid = linspace(-2, 2, 400);
y_nullcline = xgrid.^3/3 - xgrid;    % from x' = 0

% --- Plot nullclines ---
figure(1); hold on;
plot(xgrid, y_nullcline, 'k--', 'DisplayName', 'Nullcline: x''=0');

for i = 1:length(a_values)
    a = a_values(i);
    color = cmap(i, :);    
    
    % Define the ODE system as a function handle
    vdp_slowfast = @(t, X) [ ...
        X(2) - X(1)^3/3 + X(1);       % x'
        epsilon * (a - X(1))];        % y'
    
    % Integrate using ode45 or ode15s (stiff for small epsilon)
    [t, X] = ode15s(vdp_slowfast, tspan, x0);
    
    % Remove transients (first parts of simulation, where periodicity did
    % not yet occur)
    steady_idx = t > -1; % before 300
    x = X(steady_idx,1);
    y = X(steady_idx,2);
    t_steady = t(steady_idx);
    
    % --- Plot phase portrait ---
    figure(1); hold on;
    h_line = plot(x, y, 'Color', color, 'DisplayName', sprintf('a=%.2f', a));

    % ---- Mark equilibrium ----
    figure(1); hold on;
    x_eq = a;
    y_eq = x_eq^3/3 - x_eq;
    plot(x_eq, y_eq, 'o', ...
        'MarkerFaceColor', color, ...   % use line color
        'MarkerEdgeColor', color, ...   % optional, same edge
        'HandleVisibility', 'off'); 

    % --- Plot time series for x ---
    figure(2); hold on;
    plot(t_steady, x, 'DisplayName', sprintf('a=%.2f', a));

    % --- Plot time series for y ---
    figure(3); hold on;
    plot(t_steady, y, 'DisplayName', sprintf('a=%.2f', a));
    
    % --- Compute period (optional) ---
    [pks, locs] = findpeaks(x, t_steady);
    if length(locs) > 2
        period = mean(diff(locs));
        fprintf('a = %.2f --> Period ≈ %.3f\n', a, period);
    end
    
end

% Label and show plots
figure(1);
xlabel('x'); ylabel('y');
title('Phase portraits for different a');
legend show;

figure(2);
xlabel('t'); ylabel('x(t)');
title('Time series of x for different a');
legend show;

figure(3);
xlabel('t'); ylabel('y(t)');
title('Time series of y for different a');
legend show;