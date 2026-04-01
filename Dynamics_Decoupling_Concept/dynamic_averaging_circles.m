clear; clc; close all;

%% Create figure folder
figDir = 'figures_dynamics';
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

%% Parameters
N = 1000;
s = linspace(0,1,N);
eps = 0.8;
L = 2*pi;

%% Geometry
x = cos(2*pi*s);
y = sin(2*pi*s);

%% Speed profiles
v1 = ones(size(s));
v2 = 1 + eps*sin(2*pi*s);

%% Means
v_mean = 0.5*(v1 + v2);
v_harm = 2 ./ (1./v1 + 1./v2);

%% Time densities
rho1 = L ./ v1;
rho2 = L ./ v2;
rho_mean = L ./ v_mean;
rho_harm = L ./ v_harm;

%% Time maps (UNNORMALIZED)
t1 = cumtrapz(s, rho1);
t2 = cumtrapz(s, rho2);
t_mean = cumtrapz(s, rho_mean);
t_harm = cumtrapz(s, rho_harm);

%% Periods
T1 = t1(end);
T2 = t2(end);
T_mean = 0.5*(T1 + T2);
T_arith = t_mean(end);
T_harm = t_harm(end);

fprintf('\nPeriods:\n');
fprintf('T1      = %.6f\n', T1);
fprintf('T2      = %.6f\n', T2);
fprintf('T_mean  = %.6f\n', T_mean);
fprintf('T_arith = %.6f\n', T_arith);
fprintf('T_harm  = %.6f\n\n', T_harm);

%% Export helper
exportFig = @(name) exportgraphics(gcf, fullfile(figDir, name), 'ContentType','vector');

%% ---------------------------
%% Plot 1: Speed profiles
figure;
plot(s, v1, 'LineWidth',2); hold on;
plot(s, v2, 'LineWidth',2);
plot(s, v_mean, '--','LineWidth',2);
plot(s, v_harm, ':','LineWidth',2);
xlabel('s'); ylabel('speed');
title('Speed profiles');
grid on;
exportFig('speed_profiles.pdf');

%% ---------------------------
%% Plot 2: Time densities
figure;
plot(s, rho1, 'LineWidth',2); hold on;
plot(s, rho2, 'LineWidth',2);
plot(s, rho_mean, '--','LineWidth',2);
plot(s, rho_harm, ':','LineWidth',2);
xlabel('s'); ylabel('dt/ds');
title('Traversal time densities');
grid on;
exportFig('time_density.pdf');

%% ---------------------------
%% Plot 3: Time maps (UNNORMALIZED — KEY FIX)
figure;
plot(s, t1, 'LineWidth',2); hold on;
plot(s, t2, 'LineWidth',2);
plot(s, t_mean, '--','LineWidth',2);
plot(s, t_harm, ':','LineWidth',2);
xlabel('s'); ylabel('time');
title('Cumulative time maps (unnormalized)');
grid on;
exportFig('time_maps_unnormalized.pdf');

%% ---------------------------
%% OPTIONAL: normalized version (for comparison only)
t1n = t1 / T1;
t2n = t2 / T2;
tmn = t_mean / T_arith;
thn = t_harm / T_harm;

figure;
plot(s, t1n, 'LineWidth',2); hold on;
plot(s, t2n, 'LineWidth',2);
plot(s, tmn, '--','LineWidth',2);
plot(s, thn, ':','LineWidth',2);
xlabel('s'); ylabel('normalized time');
title('Cumulative time maps (normalized)');
grid on;
exportFig('time_maps_normalized.pdf');