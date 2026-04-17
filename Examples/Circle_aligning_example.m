function plot_phase_alignment_advanced
close all; clc;

set(groot, ...
    'DefaultTextInterpreter','latex', ...
    'DefaultAxesTickLabelInterpreter','latex', ...
    'DefaultAxesFontName','Times', ...
    'DefaultAxesFontSize', 12);

%% PARAMETERS
M = 400;
s = linspace(0,1,M);

% Curves
circle = [cos(2*pi*s)', sin(2*pi*s)'];
oval   = [1.5*cos(2*pi*s)', 0.7*sin(2*pi*s)'];

% Shift
tau = 0.25;
s_shift = mod(s + tau, 1);
oval_shifted = interp1(s, oval, s_shift, 'pchip');

%% ALIGNMENT
tauGrid = linspace(0,1,200);
bestErr = inf;

for t = tauGrid
    s_try = mod(s + t, 1);
    oval_try = interp1(s, oval, s_try, 'pchip');

    err = mean(sum((circle - oval_try).^2,2));
    if err < bestErr
        bestErr = err;
        bestOval = oval_try;
        bestTau = t;
    end
end

%% MEANS
mean_unaligned = 0.5*(circle + oval_shifted);
mean_aligned   = 0.5*(circle + bestOval);

%% COLOR MAP
cmap = parula(M);

%% ==========================================================
% FIGURE
%% ==========================================================
figure('Position',[100 100 600 600],'Color','w');

%% --- (1) UNALIGNED XY ---
hold on; axis equal; grid on;

for i = 1:M-1
    plot(circle(i:i+1,1), circle(i:i+1,2), 'Color', cmap(i,:), 'LineWidth',1.5);
    plot(oval_shifted(i:i+1,1), oval_shifted(i:i+1,2), '--', 'Color', cmap(i,:), 'LineWidth',1.5);
end

% starting points
plot(circle(1,1), circle(1,2), 'ko','MarkerFaceColor','b','MarkerSize',6);
plot(oval_shifted(1,1), oval_shifted(1,2), 'bo','MarkerFaceColor','b','MarkerSize',6);

plot(mean_unaligned(:,1), mean_unaligned(:,2),'Color', [0.5 0.5 0.5],'LineWidth',2);


%% EXPORT
exportgraphics(gcf,'figures_circle_aligning/circles_unaligned.pdf','ContentType','image');

%% --- (2) ALIGNED XY ---
figure('Position',[100 100 600 600],'Color','w');
hold on; axis equal; grid on;

for i = 1:M-1
    plot(circle(i:i+1,1), circle(i:i+1,2), 'Color', cmap(i,:), 'LineWidth',1.5);
    plot(bestOval(i:i+1,1), bestOval(i:i+1,2), '--', 'Color', cmap(i,:), 'LineWidth',1.5);
end

% aligned starting points
plot(circle(1,1), circle(1,2), 'ko','MarkerFaceColor','b','MarkerSize',6);
plot(bestOval(1,1), bestOval(1,2), 'bo','MarkerFaceColor','b','MarkerSize',6);

plot(mean_aligned(:,1), mean_aligned(:,2),'Color', [0.5 0.5 0.5],'LineWidth',2);


%% EXPORT
exportgraphics(gcf,'figures_circle_aligning/circles_aligned.pdf','ContentType','image');
