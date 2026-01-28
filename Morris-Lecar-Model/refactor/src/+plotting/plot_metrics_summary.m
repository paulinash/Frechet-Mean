function h = plot_metrics_summary(metrics, meanC, plotOpts)
%PLOT_METRICS_SUMMARY Basic quantitative diagnostic plots.
%
% Returns struct h with figure handles:
%   h.figDistToMean
%   h.figPairwise
%   h.figCurvature

arguments
    metrics (1,1) struct
    meanC   (:,:) double
    plotOpts struct = struct()
end

if ~isfield(plotOpts,"export"), plotOpts.export = false; end
if ~isfield(plotOpts,"outDir"),  plotOpts.outDir  = "figures"; end
plotOpts.outDir = string(plotOpts.outDir);

if plotOpts.export && ~isfolder(plotOpts.outDir)
    mkdir(plotOpts.outDir);
end

%% 1) Histograms distance vector
h = struct();
h.figHistograms = figure('Color','w','Name','Distances to mean','Position',[200 100 800 520]);
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% a) Distances to mean histogram
ax1 = nexttile(t,1);
histogram(ax1, metrics.dist.toMean, metrics.modality.edges);
grid(ax1,'on');
xlabel(ax1,'$d(C_i, \mu)$'); ylabel(ax1,'count');
title(ax1,'Distances to Frechet mean', 'Interpreter','latex');

% show modality classification
txt = sprintf('Modality: %s (peaks=%d)', string(metrics.modality.class), metrics.modality.numPeaks);
y = ylim(ax1); x = xlim(ax1);
text(ax1, x(1) + 0.05*(x(2)-x(1)), y(1) + 0.90*(y(2)-y(1)), txt);
hold(ax1,'off');

% b) Pairwise distance histogram
ax2 = nexttile(t,2);
histogram(ax2, metrics.pairwise.d, 10);
grid(ax2,'on');
xlabel(ax2,'$d(C_i, C_j)$'); ylabel(ax2,'count');
title(ax2,'Pairwise distances');

if plotOpts.export
    exportgraphics(h.figHistograms, fullfile(plotOpts.outDir, "metrics_histogram_distances.pdf"), ...
        'ContentType','image', 'Resolution', 600);
end

%% 2) Medoid vs. Mean
fig3D_Medoid_Mean = figure('Color','w','Position',[200 100 900 700]);
hold on; grid on;
plot3([meanC(:,1);meanC(1,1)], [meanC(:,2);meanC(1,2)], [meanC(:,3);meanC(1,3)], ...
    'k', 'LineWidth', 2, 'DisplayName', 'Fréchet mean');
plot3([metrics.medoid.curve(:,1);metrics.medoid.curve(1,1)], ...
    [metrics.medoid.curve(:,2);metrics.medoid.curve(1,2)], ...
    [metrics.medoid.curve(:,3);metrics.medoid.curve(1,3)], ...
    'b', 'LineWidth', 2, 'DisplayName', 'Fréchet medoid');
xlabel('V'); ylabel('w'); zlabel('y');
grid on;
title('Fréchet medoid vs. Fréchet mean');
legend();
view([-30 20]);

if plotOpts.export
    exportgraphics(fig3D_Medoid_Mean, fullfile(plotOpts.outDir, "medoid_vs_mean.pdf"), ...
        'ContentType','image', 'Resolution', 600);
end


%% 3) Curvature metrics
%  a) Plot kappa vectors for mean curve and samples curves
kappa_meanC = metrics.curvature.corr.mean.kappa;   % (M-2) x 1
S = metrics.curvature.corr.samples;                % struct array N x 1
N = numel(S);

h.figCurvature = figure();
for i = 1:N
    %subplot(N+1,1,i);
    plot(S(i).kappa, 'k-'); hold on;  % default color; can look busy
end
%subplot(N+1,1,N+1);
plot(kappa_meanC, 'r-');

grid on;
xlabel('index along curve (approx arc-length samples)');
ylabel('\kappa');
title('Curvature \kappa(s): samples (thin) and Fréchet mean (thick)');
if plotOpts.export
    exportgraphics(h.figCurvature, fullfile(plotOpts.outDir, "curvature_plots.pdf"), ...
        'ContentType','image', 'Resolution', 600);
end

%  b) Curvature metrics comparison
kappa_mean_meanC = metrics.curvature.arclen.mean.kappa_mean;
kappa_mean_samples = arrayfun(@(s) s.kappa_mean, metrics.curvature.arclen.samples);
kappa_rms_meanC = metrics.curvature.arclen.mean.kappa_rms;
kappa_rms_samples = arrayfun(@(s) s.kappa_rms, metrics.curvature.arclen.samples);
total_curvature_meanC = metrics.curvature.arclen.mean.total_curvature;
total_curvature_samples = arrayfun(@(s) s.total_curvature, metrics.curvature.arclen.samples);
bending_energy_meanC = metrics.curvature.arclen.mean.bending_energy;
bending_energy_samples = arrayfun(@(s) s.bending_energy, metrics.curvature.arclen.samples);
total_variation_meanC = metrics.curvature.arclen.mean.kappa_total_variation;
total_variation_samples = arrayfun(@(s) s.kappa_total_variation, metrics.curvature.arclen.samples);
second_difference_energy_meanC = metrics.curvature.arclen.mean.kappa_second_difference_energy;
second_difference_energy_samples = arrayfun(@(s) s.kappa_second_difference_energy, metrics.curvature.arclen.samples);

samples_list = {kappa_mean_samples, kappa_rms_samples, total_curvature_samples, bending_energy_samples, total_variation_samples, second_difference_energy_samples};
meanC_list = [kappa_mean_meanC, kappa_rms_meanC, total_curvature_meanC, bending_energy_meanC, total_variation_meanC, second_difference_energy_meanC];
title_list = {'Mean curvature', 'RMS curvature', 'Total curvature', 'Bending energy', ' total variation', 'second diff energy'};
y_label_list = {'Mean ($\kappa$)', 'RMS ($\kappa$)', '$\int \kappa ds$', '$\int \kappa^2 ds$', ' $\sum|k_i - k_{i+1}|$', '...'};

h.figCurvature_metrics = figure('Color','w','Name','Curvature summary','Position',[280 140 920 520]);
for counter=1:numel(samples_list)
    plot_curvature_samples(meanC_list(counter), samples_list{counter}, title_list{counter}, y_label_list{counter}, counter)
end

if plotOpts.export
    exportgraphics(h.figCurvature_metrics, fullfile(plotOpts.outDir, "metrics_curvature.pdf"), ...
        'ContentType','image', 'Resolution', 600);
end

end


%% Auxiliary functions
function plot_curvature_boxplot(meanC, samples, title_bp, y_label, counter)
subplot(2,3,counter);
boxplot(samples, 'Symbol', 'k.');
hold on;
plot(1, meanC, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
grid on;
title(title_bp, 'Interpreter','none');
ylabel(y_label, 'Interpreter','latex');
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax,'XTickLabel',{'samples'});
ymin = min([samples(:); meanC]);
ymax = max([samples(:); meanC]);
ylim([0.9*ymin, 1.1*ymax]);
end



function plot_curvature_samples(meanC, samples, title_bp, y_label, counter)
subplot(2,3,counter);
hold on;

% x positions (with jitter for visibility)
x = ones(size(samples));
jitter = 0.05 * randn(size(samples));
scatter(x + jitter, samples, 25, 'k', 'filled');

% plot mean in red
plot(1, meanC, 'r+', 'MarkerSize', 12, 'LineWidth', 2);

grid on;
title(title_bp, 'Interpreter','none');
ylabel(y_label, 'Interpreter','latex');

ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax,'XTick',1,'XTickLabel',{'samples'});

ymin = min([samples(:); meanC]);
ymax = max([samples(:); meanC]);
ylim([0.9*ymin, 1.1*ymax]);
xlim([0.7 1.3]);

end
