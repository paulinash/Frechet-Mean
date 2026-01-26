function print_metrics_summary(metrics)
%PRINT_METRICS_SUMMARY Pretty console summary of key diagnostics.

fprintf('\n=== Fréchet mean metrics (N=%d) ===\n', metrics.N);


% Frechet Variance vs. PairWiseDistance
fprintf('Fréchet variance: %.6g\n', metrics.dist.frechetVar);
fprintf('Average squared pairwise distance : %.6g\n', metrics.pairwise.mean);
fprintf('Ratio FV/ASPD = %.6g\n', metrics.pairwise.ratio_meanToMeanPairwise);

% Distance
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
    kappa_mean_mean = metrics.curvature.mean.kappa_mean;
    kappa_mean_samples = arrayfun(@(s) s.kappa_mean, metrics.curvature.samples);
    kappa_rms_mean = metrics.curvature.mean.kappa_rms;
    kappa_rms_samples = arrayfun(@(s) s.kappa_rms, metrics.curvature.samples);
    total_curvature_mean = metrics.curvature.mean.total_curvature;
    total_curvature_samples = arrayfun(@(s) s.total_curvature, metrics.curvature.samples);
    bending_energy_mean = metrics.curvature.mean.bending_energy;
    bending_energy_samples = arrayfun(@(s) s.bending_energy, metrics.curvature.samples);
    second_differences_mean = metrics.curvature.mean.second_diff_rms;
    second_differences_samples = arrayfun(@(s) s.second_diff_rms, metrics.curvature.samples);

    % returns mean curvature (mean kappa) of the mean curve 
    % the median of the mean curvature of the sample curves
    % the interquartile range of mean curvature of sample curves
    fprintf('Mean Curvature (mean kappa): meanCurve=%.10g, samples median=%.10g (IQR=%.10g)\n', ...
        kappa_mean_mean, median(kappa_mean_samples), iqr(kappa_mean_samples));
    % same for rms curvature
    fprintf('RMS Curvature (rms kappa): meanCurve=%.10g, samples median=%.10g (IQR=%.10g)\n', ...
        kappa_rms_mean, median(kappa_rms_samples), iqr(kappa_rms_samples));
    % same for total curvature
    fprintf('Total curvature:        meanCurve=%.10g, samples median=%.6g (IQR=%.10g)\n', ...
        total_curvature_mean, median(total_curvature_samples), iqr(total_curvature_samples));
    % same for bending energy
    fprintf('Bending energy:        meanCurve=%.10g, samples median=%.6g (IQR=%.10g)\n', ...
        bending_energy_mean, median(bending_energy_samples), iqr(bending_energy_samples));
    % same for second differences
    fprintf('Second differences:        meanCurve=%.10g, samples median=%.6g (IQR=%.10g)\n', ...
        second_differences_mean, median(second_differences_samples), iqr(second_differences_samples));
end

fprintf('==================================\n\n');
end
