function print_metrics_summary(metrics)
%PRINT_METRICS_SUMMARY Pretty console summary of key diagnostics.

fprintf('\n=== Fréchet mean metrics (N=%d) ===\n', metrics.N);


% Frechet Variance vs. PairWiseDistance
fprintf('Fréchet variance: %.6g\n', metrics.dist.frechetVar);
fprintf('Average squared pairwise distance : %.6g\n', metrics.pairwise.mean);
fprintf('Ratio FV/ASPD = %.6g\n', metrics.pairwise.ratio_meanToMeanPairwise);

% distances
fprintf('Distances to mean d(C_i,mu):\n');
fprintf('  mean = %.6g, std = %.6g, min = %.6g, max = %.6g\n', ...
    mean(metrics.dist.toMean), std(metrics.dist.toMean), ...
    min(metrics.dist.toMean), max(metrics.dist.toMean));

% modality
if isfield(metrics,'modality')
    fprintf('Distance-to-mean modality: %s (peaks=%d)\n', ...
        string(metrics.modality.class), metrics.modality.numPeaks);
end

% curvature
if isfield(metrics,'curvature')
    km = metrics.curvature.mean.kappa_rms;
    ks = arrayfun(@(s) s.kappa_rms, metrics.curvature.samples);
    bm = metrics.curvature.mean.bending_energy;
    bs = arrayfun(@(s) s.bending_energy, metrics.curvature.samples);

    fprintf('Curvature (RMS kappa): meanCurve=%.6g, samples median=%.6g (IQR=%.6g)\n', ...
        km, median(ks), iqr(ks));
    fprintf('Bending energy:        meanCurve=%.6g, samples median=%.6g (IQR=%.6g)\n', ...
        bm, median(bs), iqr(bs));
end

fprintf('==================================\n\n');
end
