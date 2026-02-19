%% Fast-Slow Curves Example: SRVF vs Arc-Length Mean
clc; clear; close all;

%% Create two synthetic fast-slow curves
num_points = 500;

% Time for curve 1
t1 = linspace(0,2*pi,num_points);
x1 = sin(t1);
y1 = cos(t1);

% Curve 2: fast first half, slow second half
t2_fast = linspace(0,pi,200).^0.5;   % fast first half (compressed time)
t2_slow = linspace(pi,2*pi,300);     % slow second half
t2 = [t2_fast, t2_slow];
x2 = sin(t2);
y2 = cos(t2);

curves = { [x1; y1], [x2; y2] };
times  = { t1, t2 };

%% ---------- 1) Arc-Length Mean (pure shape) ----------
arc_resampled = zeros(2, num_points, 2);

for i = 1:2
    c = curves{i}; % 2 × N
    
    % Compute arc-length parameter s(t)
    ds = sqrt(sum(diff(c,1,2).^2,1));   % 1 × (N-1)
    s = [0, cumsum(ds)];                % 1 × N
    s = s / s(end);                     % normalize to [0,1]
    
    % Resample x(s), y(s) onto a uniform arc-length grid
    s_uniform = linspace(0,1,num_points);
    
    arc_resampled(:, :, i) = interp1(s, c', s_uniform, 'linear')';
end

% Arc-length mean shape
arc_mean_shape = mean(arc_resampled, 3);  % 2 × num_points


%% ---------- 2) SRVF Mean (elastic time warping) ----------
% Initialize SRVF-mean using curve 1
mean_x = x1;
mean_y = y1;

num_iter = 5;  % small number is enough for 2 curves

for iter = 1:num_iter
    aligned = zeros(2, num_points, 2);
    
    for i = 1:2
        c = curves{i}; % 2×N
        
        % arc-length of curve
        ds_c = sqrt(sum(diff(c,1,2).^2,1));
        s_c = [0, cumsum(ds_c)];
        s_c = s_c / s_c(end);
        
        % arc-length of current template mean
        ds_mean = sqrt(diff(mean_x).^2 + diff(mean_y).^2);
        s_mean = [0, cumsum(ds_mean)];
        s_mean = s_mean / s_mean(end);
        
        % ensure uniqueness for interpolation
        [s_unique, ia] = unique(s_c);

        % warp curve to mean-parameter domain
        aligned(:,:,i) = interp1(s_unique, c(:,ia)', s_mean, 'linear', 'extrap')';
    end
    
    % Update mean
    mean_x = mean(aligned(1,:,:),3);
    mean_y = mean(aligned(2,:,:),3);
end

srvf_mean_shape = [mean_x; mean_y];

%% ---------- 3) SRVF mean resampled back to original t ----------
% Resample the SRVF mean along the t1 grid
t_uniform = linspace(0,1,num_points);

ds_mean = sqrt(diff(mean_x).^2 + diff(mean_y).^2);
s_mean = [0, cumsum(ds_mean)];
s_mean = s_mean / s_mean(end);

% Resample according to *original time grid*
mean_x_resampled = interp1(s_mean, mean_x, t_uniform, 'linear', 'extrap');
mean_y_resampled = interp1(s_mean, mean_y, t_uniform, 'linear', 'extrap');


%% ---------- Plotting ----------
figure; hold on; axis equal; grid on;
title('SRVF Mean vs Arc-Length Mean vs Fast-Slow Dynamics');
xlabel('x'); ylabel('y');

plot(x1,y1,'b','LineWidth',1.5);
plot(x2,y2,'r','LineWidth',1.5);

plot(arc_mean_shape(1,:), arc_mean_shape(2,:), 'g', 'LineWidth', 3);
plot(srvf_mean_shape(1,:), srvf_mean_shape(2,:), 'k--', 'LineWidth', 3);
plot(mean_x_resampled, mean_y_resampled, 'm', 'LineWidth', 3);

legend('Curve 1','Curve 2',...
       'Arc-length mean (no timing)', ...
       'SRVF elastic mean (warped timing)', ...
       'SRVF mean w/ original timing preserved');

%% ---------- Animation Setup ----------
figure; hold on; axis equal; grid on;
xlabel('x'); ylabel('y');
title('Animated Trajectories: Fast-Slow, Arc Mean, SRVF Mean');

% Plot full curves faintly in background
plot(x1,y1,'b:');
plot(x2,y2,'r:');
plot(arc_mean_shape(1,:), arc_mean_shape(2,:), 'g:');
plot(srvf_mean_shape(1,:), srvf_mean_shape(2,:), 'k:');
plot(mean_x_resampled, mean_y_resampled, 'm:');

legend('Curve 1 (static)','Curve 2 (static)',...
       'Arc-length mean (static)',...
       'SRVF elastic mean (static)',...
       'SRVF mean w/ timing (static)',...
       'Location','BestOutside');

% Create animated markers
h1 = plot(x1(1), y1(1), 'bo', 'MarkerSize', 10, 'MarkerFaceColor','b');
h2 = plot(x2(1), y2(1), 'ro', 'MarkerSize', 10, 'MarkerFaceColor','r');
h_arc  = plot(arc_mean_shape(1,1),  arc_mean_shape(2,1),  'go', 'MarkerSize', 10, 'MarkerFaceColor','g');
h_srvf = plot(srvf_mean_shape(1,1), srvf_mean_shape(2,1), 'ko', 'MarkerSize', 10, 'MarkerFaceColor','k');
h_srvf_t = plot(mean_x_resampled(1), mean_y_resampled(1),'mo','MarkerSize',10,'MarkerFaceColor','m');

drawnow;

%% ---------- Animate ----------
for k = 1:num_points
    
    % Move each marker according to its own time parameterization
    set(h1, 'XData', x1(k), 'YData', y1(k));
    set(h2, 'XData', x2(k), 'YData', y2(k));
    set(h_arc, 'XData', arc_mean_shape(1,k), 'YData', arc_mean_shape(2,k));
    set(h_srvf, 'XData', srvf_mean_shape(1,k), 'YData', srvf_mean_shape(2,k));
    set(h_srvf_t, 'XData', mean_x_resampled(k), 'YData', mean_y_resampled(k));

    drawnow;
    pause(0.01);  % Adjust playback speed
end