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

% curvature arrays
ks = arrayfun(@(s) s.kappa_rms, metrics.curvature.samples);
bs = arrayfun(@(s) s.bending_energy, metrics.curvature.samples);

T = table();
T.N = metrics.N;

T.d_to_mean_mean = mean(d);
T.d_to_mean_std  = std(d);
T.d_to_mean_min  = min(d);
T.d_to_mean_max  = max(d);

T.frechetVar = metrics.dist.frechetVar;
T.pairwise_mean = metrics.pairwise.mean;
T.ratio_meanToMeanPairwise = metrics.pairwise.ratio_meanToMeanPairwise;


T.modality_class = string(metrics.modality.class);
T.modality_peaks = metrics.modality.numPeaks;

T.kappa_rms_meanCurve = metrics.curvature.mean.kappa_rms;
T.kappa_rms_medianSamples = median(ks);
T.kappa_rms_iqrSamples = iqr(ks);

T.bend_meanCurve = metrics.curvature.mean.bending_energy;
T.bend_medianSamples = median(bs);
T.bend_iqrSamples = iqr(bs);



if ~plotOpts.export
    return;
end

% --- export CSV ---
writetable(T, fullfile(plotOpts.outDir, "metrics_summary.csv"));

end
