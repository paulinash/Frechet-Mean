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


T = table();
T.N = metrics.N;

T.d_to_mean_mean = mean(d);
T.d_to_mean_std  = std(d);
T.d_to_mean_min  = min(d);
T.d_to_mean_max  = max(d);

T.frechetVar = metrics.dist.frechetVar;
T.pairwise_mean = metrics.pairwise.mean;
T.ratio_meanToMeanPairwise = metrics.pairwise.ratio_meanToMeanPairwise;

%% Medoid
T.frechetVar_medoid = metrics.medoid.cost;
T.distance_Fmedoid_Fmean = metrics.medoid.d_to_mean;
T.ratio_FVmedoid_FVmean = metrics.medoid.varRatio;

%% Modality
T.modality_class = string(metrics.modality.class);
T.modality_peaks = metrics.modality.numPeaks;

%% Curvature
%curvature stats of mean curve
T.kappa_mean_meanC = metrics.curvature.arclen.mean.kappa_mean;
T.kappa_rms_meanC = metrics.curvature.arclen.mean.kappa_rms;
T.total_curvature_meanC = metrics.curvature.arclen.mean.total_curvature;
T.bending_energy_meanC = metrics.curvature.arclen.mean.bending_energy;
T.total_variation_meanC = metrics.curvature.arclen.mean.kappa_total_variation;
T.second_difference_energy_meanC = metrics.curvature.arclen.mean.kappa_second_difference_energy;

% --- curvature sample curves arrays ---
kappa_mean_samples = arrayfun(@(s) s.kappa_mean, metrics.curvature.arclen.samples);
kappa_rms_samples = arrayfun(@(s) s.kappa_rms, metrics.curvature.arclen.samples);
total_curvature_samples = arrayfun(@(s) s.total_curvature, metrics.curvature.arclen.samples);
bending_energy_samples = arrayfun(@(s) s.bending_energy, metrics.curvature.arclen.samples);
total_variation_samples = arrayfun(@(s) s.kappa_total_variation, metrics.curvature.arclen.samples);
second_diff_energy_samples = arrayfun(@(s) s.kappa_second_difference_energy, metrics.curvature.arclen.samples);

% mean of curvature stats of sample curves
T.kappa_mean_meanSamples = mean(kappa_mean_samples);
T.kappa_rms_meanSamples = mean(kappa_rms_samples);
T.total_curvature_meanSamples = mean(total_curvature_samples);
T.bending_energy_meanSamples = mean(bending_energy_samples);
T.total_variation_meanSamples = mean(total_variation_samples);
T.second_difference_energy_meanSamples = mean(second_diff_energy_samples);

% Ratios
T.kappa_mean_ratio = T.kappa_mean_meanC/T.kappa_mean_meanSamples;
T.kappa_rms_ratio = T.kappa_rms_meanC / T.kappa_rms_meanSamples;
T.total_curvature_ratio = T.total_curvature_meanC  / T.total_curvature_meanSamples;
T.bending_energy_ratio = T.bending_energy_meanC / T.bending_energy_meanSamples;
T.total_variation_ratio = T.total_variation_meanC / T.total_variation_meanSamples;
T.second_diff_ratio = T.second_difference_energy_meanC / T.second_difference_energy_meanSamples;



if ~plotOpts.export
    return;
end

% --- export CSV ---
writetable(T, fullfile(plotOpts.outDir, "metrics_summary.csv"));

end
