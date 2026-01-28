function Tlong = export_metrics_table(metrics, plotOpts)
%EXPORT_METRICS_TABLE Export key metrics summary as a grouped Excel table.
%
% Output:
%   - writes metrics_summary.xlsx in plotOpts.outDir (default "figures")
% Returns:
%   Tlong: table with columns Group, Name, Value

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

% --- row builders ---
Group = strings(0,1);
Name  = strings(0,1);
Value = strings(0,1);

add = @(g,n,v) addRow(g,n,v);

    function addRow(g,n,v)
        Group(end+1,1) = string(g);
        Name(end+1,1)  = string(n);
        Value(end+1,1) = toStr(v);
    end

    function s = toStr(v)
        if isstring(v) || ischar(v)
            s = string(v);
        elseif islogical(v)
            s = string(v);
        elseif isnumeric(v)
            if isscalar(v)
                s = string(v);
            else
                % keep non-scalars compact
                s = string(mat2str(v, 6));
            end
        else
            s = string(v);
        end
    end

addHeader = @(g) addRow(g, upper(string(g)), "");  % group header row
addBlank  = @() addRow("", "", "");                % blank separator row

% --- build rows ---
d = metrics.dist.toMean;

% BASIC
addHeader("Basic");
add("Basic","N", metrics.N);
add("Basic","d_to_mean_mean", mean(d));
add("Basic","d_to_mean_std",  std(d));
add("Basic","d_to_mean_min",  min(d));
add("Basic","d_to_mean_max",  max(d));
add("Basic","frechetVar", metrics.dist.frechetVar);
add("Basic","pairwise_mean_sq", metrics.pairwise.mean);
add("Basic","ratio_meanToMeanPairwise", metrics.pairwise.ratio_meanToMeanPairwise);
addBlank();

% MEDOID
addHeader("Medoid");
add("Medoid","medoid_idx", metrics.medoid.idx);
add("Medoid","frechetVar_medoid", metrics.medoid.cost);
add("Medoid","distance_Fmedoid_Fmean", metrics.medoid.d_to_mean);
add("Medoid","ratio_FVmedoid_FVmean", metrics.medoid.varRatio);
addBlank();

% MODALITY
addHeader("Modality");
add("Modality","modality_class", string(metrics.modality.class));
add("Modality","modality_peaks", metrics.modality.numPeaks);
addBlank();

% CURVATURE
addHeader("Curvature");
add("Curvature","kappa_mean_meanC", metrics.curvature.arclen.mean.kappa_mean);
add("Curvature","kappa_rms_meanC", metrics.curvature.arclen.mean.kappa_rms);
add("Curvature","total_curvature_meanC", metrics.curvature.arclen.mean.total_curvature);
add("Curvature","bending_energy_meanC", metrics.curvature.arclen.mean.bending_energy);
add("Curvature","total_variation_meanC", metrics.curvature.arclen.mean.kappa_total_variation);
add("Curvature","second_difference_energy_meanC", metrics.curvature.arclen.mean.kappa_second_difference_energy);

% sample stats means
S = metrics.curvature.arclen.samples;
kappa_mean_samples = arrayfun(@(s) s.kappa_mean, S);
kappa_rms_samples  = arrayfun(@(s) s.kappa_rms,  S);
total_curv_samples = arrayfun(@(s) s.total_curvature, S);
bendE_samples      = arrayfun(@(s) s.bending_energy, S);
totVar_samples     = arrayfun(@(s) s.kappa_total_variation, S);
secDiff_samples    = arrayfun(@(s) s.kappa_second_difference_energy, S);

kmS = mean(kappa_mean_samples);
krS = mean(kappa_rms_samples);
tcS = mean(total_curv_samples);
beS = mean(bendE_samples);
tvS = mean(totVar_samples);
sdS = mean(secDiff_samples);

add("Curvature","kappa_mean_meanSamples", kmS);
add("Curvature","kappa_rms_meanSamples",  krS);
add("Curvature","total_curvature_meanSamples", tcS);
add("Curvature","bending_energy_meanSamples",  beS);
add("Curvature","total_variation_meanSamples", tvS);
add("Curvature","second_difference_energy_meanSamples", sdS);

% ratios
add("Curvature","kappa_mean_ratio", metrics.curvature.arclen.mean.kappa_mean / kmS);
add("Curvature","kappa_rms_ratio",  metrics.curvature.arclen.mean.kappa_rms  / krS);
add("Curvature","total_curvature_ratio", metrics.curvature.arclen.mean.total_curvature / tcS);
add("Curvature","bending_energy_ratio",  metrics.curvature.arclen.mean.bending_energy / beS);
add("Curvature","total_variation_ratio", metrics.curvature.arclen.mean.kappa_total_variation / tvS);
add("Curvature","second_diff_ratio",     metrics.curvature.arclen.mean.kappa_second_difference_energy / sdS);

% --- output table ---
Tlong = table(Group, Name, Value);

if ~plotOpts.export
    return;
end

outPath = fullfile(plotOpts.outDir, "metrics_summary.xlsx");

% Write to Excel
writetable(Tlong, outPath, 'FileType','spreadsheet', 'Sheet', 'Summary');

end
