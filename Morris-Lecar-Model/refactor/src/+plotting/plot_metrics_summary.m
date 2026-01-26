function h = plot_metrics_summary(metrics, plotOpts)
%PLOT_METRICS_SUMMARY Basic quantitative diagnostic plots.
%
% Returns struct h with figure handles:
%   h.figDistToMean
%   h.figPairwise
%   h.figCurvature

arguments
    metrics (1,1) struct
    plotOpts struct = struct()
end

if ~isfield(plotOpts,"export"), plotOpts.export = false; end
if ~isfield(plotOpts,"outDir"),  plotOpts.outDir  = "figures"; end
plotOpts.outDir = string(plotOpts.outDir);

if plotOpts.export && ~isfolder(plotOpts.outDir)
    mkdir(plotOpts.outDir);
end

h = struct();

%% 1) Distances to mean histogram
h.figDistToMean = figure('Color','w','Name','Distances to mean','Position',[200 100 800 520]);
histogram(metrics.dist.toMean, metrics.modality.edges);
grid on;
xlabel('$d(C_i, \mu)$'); ylabel('count');
title('Distances to Frechet mean', 'Interpreter','latex');

% show modality classification
txt = sprintf('Modality: %s (peaks=%d)', string(metrics.modality.class), metrics.modality.numPeaks);
y = ylim; x = xlim;
text(x(1) + 0.05*(x(2)-x(1)), y(1) + 0.90*(y(2)-y(1)), txt);

if plotOpts.export
    exportgraphics(h.figDistToMean, fullfile(plotOpts.outDir, "metrics_dist_to_mean.pdf"), ...
        'ContentType','image', 'Resolution', 600);
end


%% 2) Pairwise distance histogram
h.figPairwise = figure('Color','w','Name','Pairwise distances','Position',[240 120 800 520]);
histogram(metrics.pairwise.d, 30);
grid on;
xlabel('$d(C_i, C_j)$'); ylabel('count');
title('Pairwise distances');

if isfield(metrics.pairwise,'opts') && isfield(metrics.pairwise.opts,'allowShift')
    if metrics.pairwise.opts.allowShift
        subtitle('shift-invariant (min over \tau)');
    else
        subtitle('as-aligned (no shift search)');
    end
end

if plotOpts.export
    exportgraphics(h.figPairwise, fullfile(plotOpts.outDir, "metrics_pairwise_dist.pdf"), ...
        'ContentType','image', 'Resolution', 600);
end


%% plot kappa
% Plot kappa vectors: all samples + mean
kappa_meanC = metrics.curvature.mean.kappa;   % (M-2) x 1
S = metrics.curvature.samples;                % struct array N x 1
N = numel(S);

figure();
% plot sample kappas
for i = 1:N
    subplot(N,1,i);
    plot(S(i).kappa, 'k-');   % default color; can look busy
end

% plot mean kappa on top
subplot(N,1,N);
plot(kappa_meanC, 'r-', 'LineWidth', 2);

grid on;
xlabel('index along curve (approx arc-length samples)');
ylabel('\kappa');
title('Curvature \kappa(s): samples (thin) and Fréchet mean (thick)');


%% 3) Curvature metrics comparison
N = metrics.N;
kappa_mean_meanC = metrics.curvature.mean.kappa_mean;
kappa_mean_samples = arrayfun(@(s) s.kappa_mean, metrics.curvature.samples);
kappa_rms_meanC = metrics.curvature.mean.kappa_rms;
kappa_rms_samples = arrayfun(@(s) s.kappa_rms, metrics.curvature.samples);
total_curvature_meanC = metrics.curvature.mean.total_curvature;
total_curvature_samples = arrayfun(@(s) s.total_curvature, metrics.curvature.samples);
bending_energy_meanC = metrics.curvature.mean.bending_energy;
bending_energy_samples = arrayfun(@(s) s.bending_energy, metrics.curvature.samples);
second_differences_meanC = metrics.curvature.mean.second_diff_rms;
second_differences_samples = arrayfun(@(s) s.second_diff_rms, metrics.curvature.samples);

samples_list = {kappa_mean_samples, kappa_rms_samples, total_curvature_samples, bending_energy_samples, second_differences_samples};
meanC_list = [kappa_mean_meanC, kappa_rms_meanC, total_curvature_meanC, bending_energy_meanC, second_differences_meanC];
title_list = {'Mean curvature', 'RMS curvature', 'Total curvature', 'Bending energy', 'Second differences'};
y_label_list = {'Mean ($\kappa$)', 'RMS ($\kappa$)', '$\int \kappa ds$', '$\int \kappa^2 ds$', 'Second differences RMS'};

h.figCurvature = figure('Color','w','Name','Curvature summary','Position',[280 140 920 520]);
for counter=1:numel(samples_list)
    plot_curvature_boxplot(meanC_list(counter), samples_list{counter}, title_list{counter}, y_label_list{counter}, counter)
end

if plotOpts.export
    exportgraphics(h.figCurvature, fullfile(plotOpts.outDir, "metrics_curvature.pdf"), ...
        'ContentType','image', 'Resolution', 600);
end

end

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
