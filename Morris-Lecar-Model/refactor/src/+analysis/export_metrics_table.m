function T = export_metrics_table(metrics, plotOpts)
%EXPORT_METRICS_TABLE Export key metrics summary as CSV and a table figure (PDF/PNG).
%
% Returns table T.

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

% --- build a compact summary table (1 row) ---
d = metrics.dist.toMean;
pw = metrics.pairwise.d;

% curvature sample curves arrays
kappa_mean_samples = arrayfun(@(s) s.kappa_mean, metrics.curvature.samples);
kappa_rms_samples = arrayfun(@(s) s.kappa_rms, metrics.curvature.samples);
total_curvature_samples = arrayfun(@(s) s.total_curvature, metrics.curvature.samples);
bending_energy_samples = arrayfun(@(s) s.bending_energy, metrics.curvature.samples);
second_differences_samples = arrayfun(@(s) s.second_diff_rms, metrics.curvature.samples);


T = table();
T.N = metrics.N;

T.d_to_mean_mean = mean(d);
T.d_to_mean_std  = std(d);
T.d_to_mean_min  = min(d);
T.d_to_mean_max  = max(d);

T.frechetVar = metrics.dist.frechetVar;
T.pairwise_mean = metrics.pairwise.mean;
T.ratio_meanToMeanPairwise = metrics.pairwise.ratio_meanToMeanPairwise;

% modality
T.modality_class = string(metrics.modality.class);
T.modality_peaks = metrics.modality.numPeaks;

%curvature stats of mean curve
T.kappa_mean_meanC = metrics.curvature.mean.kappa_mean;
T.kappa_rms_meanC = metrics.curvature.mean.kappa_rms;
T.total_curvature_meanC = metrics.curvature.mean.total_curvature;
T.bending_energy_meanC = metrics.curvature.mean.bending_energy;
T.second_diff_meanC = metrics.curvature.mean.second_diff_rms;

% median of curvature stats of sample curves
T.kappa_mean_medianSamples = median(kappa_mean_samples);
T.kappa_rms_medianSamples = median(kappa_rms_samples);
T.total_curvature_medianSamples = median(total_curvature_samples);
T.bending_energy_medianSamples = median(bending_energy_samples);
T.second_differences_medianSamples = median(second_differences_samples);

% IQR of curvature stats of samples curves
T.kappa_mean_iqrSamples = iqr(kappa_mean_samples);
T.kappa_rms_iqrSamples = iqr(kappa_rms_samples);
T.total_curvature_iqrSamples = iqr(total_curvature_samples);
T.bending_energy_iqrSamples = iqr(bending_energy_samples);
T.second_differences_iqrSamples = iqr(second_differences_samples);


if ~plotOpts.export
    return;
end

% --- export CSV ---
writetable(T, fullfile(plotOpts.outDir, "metrics_summary.csv"));

end
