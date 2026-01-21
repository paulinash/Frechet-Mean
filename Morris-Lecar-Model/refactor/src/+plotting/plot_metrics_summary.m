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
title('Distances to Fréchet mean', 'Interpreter','latex');

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


%% 3) Curvature comparison
N = metrics.N;
kMean = metrics.curvature.mean.kappa_rms;

kSamp = zeros(N,1);
bendSamp = zeros(N,1);
for i = 1:N
    kSamp(i) = metrics.curvature.samples(i).kappa_rms;
    bendSamp(i) = metrics.curvature.samples(i).bending_energy;
end

h.figCurvature = figure('Color','w','Name','Curvature summary','Position',[280 140 920 520]);

subplot(1,2,1);
boxplot(kSamp);
hold on;
plot(1, kMean, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
grid on;
title('RMS curvature');
ylabel('RMS($\kappa$)');
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax,'XTickLabel',{'samples'});

subplot(1,2,2);
boxplot(bendSamp);
hold on;
plot(1, metrics.curvature.mean.bending_energy, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
grid on;
title('Bending energy');
ylabel('$\int \kappa^2 ds$');
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax,'XTickLabel',{'samples'});

if plotOpts.export
    exportgraphics(h.figCurvature, fullfile(plotOpts.outDir, "metrics_curvature.pdf"), ...
        'ContentType','image', 'Resolution', 600);
end

end
