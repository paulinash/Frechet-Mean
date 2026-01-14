%% Van der Pol System - SRVF Mean with Zero-Crossing Period Detection

clc; clear; close all;

%% Parameters
epsilon = 0.1;
a_values = [-0.8, 0.8];      % Example a values, including a=0
T_sim = 200;                     % Total simulation time
dt = 0.01;                       
transient_time = 100;             % Time to reach limit cycle
num_points = 500;                 % Points per period for resampling

num_a = length(a_values);
all_x = cell(num_a,1);
all_y = cell(num_a,1);
all_t = cell(num_a,1);

%% Simulation loop
for i = 1:num_a
    a = a_values(i);
    vdp = @(t,Y) [Y(2) - Y(1)^3/3 + Y(1); epsilon*(a - Y(1))];
    Y0 = [0;0];
    tspan = 0:dt:T_sim;
    [t,Y] = ode45(vdp, tspan, Y0);
    
    % Remove transient
    idx = t >= transient_time;
    t_lc = t(idx);
    x_lc = Y(idx,1);
    y_lc = Y(idx,2);
    
    % --- Detect first full period using zero-crossings ---
    cross_idx = find(x_lc(1:end-1) < 0 & x_lc(2:end) >= 0);  % zero-crossings upward
    if length(cross_idx) < 2
        error('Not enough zero-crossings to determine a period for a = %.2f', a);
    end
    
    t_period = t_lc(cross_idx(1):cross_idx(2));
    x_period = x_lc(cross_idx(1):cross_idx(2));
    y_period = y_lc(cross_idx(1):cross_idx(2));
    
    % Resample to uniform grid along its own period
    t_uniform = linspace(t_period(1), t_period(end), num_points);
    x_uniform = interp1(t_period, x_period, t_uniform, 'spline');
    y_uniform = interp1(t_period, y_period, t_uniform, 'spline');
    
    all_x{i} = x_uniform;
    all_y{i} = y_uniform;
    all_t{i} = t_uniform;
end

%% --- SRVF / Elastic Mean computation (shape only) ---
mean_x = all_x{1};
mean_y = all_y{1};
num_iter = 5;

for iter = 1:num_iter
    aligned_curves = zeros(num_points, num_a, 2); % x and y
    for i = 1:num_a
        c_x = all_x{i};
        c_y = all_y{i};
        
        % Arc-length parameterization of curve i
        ds = sqrt(diff(c_x).^2 + diff(c_y).^2);
        s_curve = [0, cumsum(ds)];
        s_curve = s_curve / s_curve(end);
        
        % Arc-length of current mean
        ds_mean = sqrt(diff(mean_x).^2 + diff(mean_y).^2);
        s_mean = [0, cumsum(ds_mean)];
        s_mean = s_mean / s_mean(end);
        
        % Ensure unique points for interpolation
        [s_unique, ia] = unique(s_curve);
        x_unique = c_x(ia);
        y_unique = c_y(ia);
        
        % Interpolate onto mean arc-length
        aligned_curves(:,i,1) = interp1(s_unique, x_unique, s_mean, 'linear', 'extrap');
        aligned_curves(:,i,2) = interp1(s_unique, y_unique, s_mean, 'linear', 'extrap');
    end
    
    % Update mean curve in shape space
    mean_x = mean(aligned_curves(:,:,1),2)';
    mean_y = mean(aligned_curves(:,:,2),2)';
end

%% --- Resample mean along original time grid of first curve ---
t_mean = linspace(all_t{1}(1), all_t{1}(end), num_points);
ds_mean = sqrt(diff(mean_x).^2 + diff(mean_y).^2);
s_mean = [0, cumsum(ds_mean)];
s_mean = s_mean / s_mean(end);
t_grid_norm = linspace(0,1,num_points);
mean_x_resampled = interp1(s_mean, mean_x, t_grid_norm, 'linear', 'extrap');
mean_y_resampled = interp1(s_mean, mean_y, t_grid_norm, 'linear', 'extrap');

%% --- Plot phase space ---
figure; hold on;
colors = lines(num_a);
for i = 1:num_a
    plot(all_x{i}, all_y{i}, 'Color', colors(i,:), 'LineWidth',1.5);
end
plot(mean_x_resampled, mean_y_resampled, 'k', 'LineWidth',3);
xlabel('x'); ylabel('y'); title('Van der Pol - SRVF Mean with Original Timing');
grid on; axis equal;

%% --- Plot x(t) ---
figure; hold on;
for i = 1:num_a
    plot(all_t{i}, all_x{i}, 'Color', colors(i,:), 'LineWidth',1.5);
end
plot(all_t{1}, mean_x_resampled, 'k', 'LineWidth',3);
xlabel('Time'); ylabel('x(t)'); title('x(t) trajectories and SRVF Mean'); grid on;

%% --- Plot y(t) ---
figure; hold on;
for i = 1:num_a
    plot(all_t{i}, all_y{i}, 'Color', colors(i,:), 'LineWidth',1.5);
end
plot(all_t{1}, mean_y_resampled, 'k', 'LineWidth',3);
xlabel('Time'); ylabel('y(t)'); title('y(t) trajectories and SRVF Mean'); grid on;

%% --- Optional animation along mean ---
figure; hold on;
for i = 1:num_a
    plot(all_x{i}, all_y{i}, 'Color', colors(i,:), 'LineWidth',1.5);
end
h_mean = plot(mean_x_resampled, mean_y_resampled, 'k', 'LineWidth',3);
xlabel('x'); ylabel('y'); title('Animation along SRVF Mean'); grid on; axis equal;

markers = gobjects(num_a+1,1);
for i = 1:num_a
    markers(i) = plot(all_x{i}(1), all_y{i}(1),'o','MarkerSize',8,'MarkerFaceColor',colors(i,:),'MarkerEdgeColor','k');
end
markers(end) = plot(mean_x_resampled(1), mean_y_resampled(1),'o','MarkerSize',10,'MarkerFaceColor','k','MarkerEdgeColor','k');

for k = 1:num_points
    for i = 1:num_a
        markers(i).XData = all_x{i}(k);
        markers(i).YData = all_y{i}(k);
    end
    markers(end).XData = mean_x_resampled(k);
    markers(end).YData = mean_y_resampled(k);
    drawnow; pause(0.01);
end