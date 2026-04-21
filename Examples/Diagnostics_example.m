clear; close all; clc;

set(groot, ...
    'DefaultTextInterpreter','latex', ...
    'DefaultAxesTickLabelInterpreter','latex', ...
    'DefaultAxesFontName','Times', ...
    'DefaultAxesFontSize', 12);

%% PARAMETERS
N1 = 20; % circles
N2 = 15; % double loops
N = N1 + N2;

M = 400;
s = linspace(0,1,M);

maxIter = 20;

%% GENERATE DATA

curves = zeros(N,2,M);

% --- circles ---
for i = 1:N1
    r = 1 + 0.1*randn;
    curves(i,1,:) = r*cos(2*pi*s);
    curves(i,2,:) = r*sin(2*pi*s);
end

% --- double loops (true figure-eight / double loop) ---
for i = 1:N2
    idx = N1 + i;
    
    a = 1 + 0.05*randn; % small variation
    
    curves(idx,1,:) = a * sin(2*pi*s);
    curves(idx,2,:) = a * sin(4*pi*s);
end

%% DISTANCE FUNCTION

dist_curve = @(c1,c2) min_shift_distance(c1,c2);

function d = min_shift_distance(c1,c2)
    M = size(c1,2);
    dmin = inf;
    for tau = 0:M-1
        c2_shift = circshift(c2,[0 tau]);
        dtmp = mean(sum((c1 - c2_shift).^2,1));
        dmin = min(dmin,dtmp);
    end
    d = dmin;
end

%% FRECHET MEAN

mu = squeeze(mean(curves,1));

for it = 1:maxIter
    aligned = zeros(size(curves));
    
    for k = 1:N
        best_shift = 0;
        best_val = inf;
        
        for tau = 0:M-1
            c_shift = circshift(squeeze(curves(k,:,:)),[0 tau]);
            val = mean(sum((mu - c_shift).^2,1));
            if val < best_val
                best_val = val;
                best_shift = tau;
            end
        end
        
        aligned(k,:,:) = circshift(squeeze(curves(k,:,:)),[0 best_shift]);
    end
    
    mu = squeeze(mean(aligned,1));
end

%% MEDOID

D = zeros(N);
for i=1:N
    for j=1:N
        D(i,j) = dist_curve(squeeze(curves(i,:,:)), squeeze(curves(j,:,:)));
    end

end

[~, med_idx] = min(sum(D,2));
medoid = squeeze(curves(med_idx,:,:));

%% VARIANCES + DISTANCES

VarF = mean(arrayfun(@(k) dist_curve(mu, squeeze(curves(k,:,:))),1:N));
VarMed = mean(arrayfun(@(k) dist_curve(medoid, squeeze(curves(k,:,:))),1:N));
mean_med_dist = dist_curve(mu,medoid);

delta_F = min(arrayfun(@(k) dist_curve(mu, squeeze(curves(k,:,:))),1:N));

delta_med = inf;
for k=1:N
    if k ~= med_idx
        dtmp = dist_curve(medoid, squeeze(curves(k,:,:)));
        delta_med = min(delta_med,dtmp);
    end
end

S = delta_F / delta_med;

%% CURVATURE

function [TC,Ebend] = curvature_metrics(c)
    x = c(1,:); y = c(2,:);
    
    dx = gradient(x); dy = gradient(y);
    ddx = gradient(dx); ddy = gradient(dy);
    
    kappa = abs(dx.*ddy - dy.*ddx) ./ (dx.^2 + dy.^2).^(3/2);
    kappa(isnan(kappa)) = 0;
    
    TC = trapz(kappa);
    Ebend = trapz(kappa.^2);
end

TCs = zeros(N,1);
Ebends = zeros(N,1);

for k=1:N
    [TCs(k), Ebends(k)] = curvature_metrics(squeeze(curves(k,:,:)));
end

[TC_mu, Ebend_mu] = curvature_metrics(mu);

%% =======================
%% PLOTTING FUNCTIONS
%% =======================

function plot_curvature_samples(meanC, samples, filename)
    figure('Color','w','Position',[100 100 600 600]); hold on;

    x = ones(size(samples));
    jitter = 0.01 * randn(size(samples));
    scatter(x + jitter, samples, 25, 'k', 'filled');

    plot(1, meanC, 'r+', 'MarkerSize', 12, 'LineWidth', 2);

    ax = gca;
    ax.YGrid = 'on';
    ax.XGrid = 'off';
    ax.GridAlpha = 0.15;
    ax.LineWidth = 1;
    ax.XTick = [];

    ymin = min([samples(:); meanC]);
    ymax = max([samples(:); meanC]);
    ylim([0.9*ymin, 1.1*ymax]);
    xlim([0.9 1.1]);

    exportgraphics(gcf,filename,'ContentType','image');

end

%% =======================
%% PLOTS
%% =======================
% --- curvature plots ---
plot_curvature_samples(TC_mu, TCs, 'figures_diagnostics/TC_plot.pdf');
plot_curvature_samples(Ebend_mu, Ebends, 'figures_diagnostics/bending_energy_plot.pdf');

% --- circles ---
figure('Color','w','Position',[100 100 600 600]); hold on; grid on;
axis([-1.5 1.5 -1.5 1.5])
for k=1:N1
    plot(squeeze(curves(k,1,:)), squeeze(curves(k,2,:)));
end
exportgraphics(gcf,'figures_diagnostics/circles.pdf','ContentType','image');

% --- double loops ---
figure('Color','w','Position',[100 100 600 600]); hold on; grid on;
axis([-1.5 1.5 -1.5 1.5])
for k=N1+1:N
    plot(squeeze(curves(k,1,:)), squeeze(curves(k,2,:)));
end
exportgraphics(gcf,'figures_diagnostics/double_loops.pdf','ContentType','image');


% --- medoid and mean ---
figure('Color','w','Position',[100 100 600 600]); hold on; grid on;
axis([-1.5 1.5 -1.5 1.5])
plot(mu(1,:),mu(2,:),'r','LineWidth',2);
plot(medoid(1,:),medoid(2,:),'b','LineWidth',2);
exportgraphics(gcf,'figures_diagnostics/medoid_mean.pdf','ContentType','image');

%% =======================
%% PRINT NUMBERS
%% =======================

fprintf('--- Diagnostics ---\n');
fprintf('Var_F = %.4f\n',VarF);
fprintf('Var_med = %.4f\n',VarMed);
fprintf('Mean-Medoid Dist = %.4f\n',mean_med_dist);
fprintf('delta_F = %.4f\n',delta_F);
fprintf('delta_med = %.7f\n',delta_med);
fprintf('Isolation S = %.4f\n',S);