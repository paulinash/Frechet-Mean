function out = exp_ml_frechet_mean()
%EXP_ML_FRECHET_MEAN Morris–Lecar bursting + (shift-aligned) Fréchet mean.
%
%
close all; clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(projectRoot,'src')));

% Get parameters
S = config_exp_ml_frechet_mean(string(projectRoot));

%% Model parameters
p0 = ml.model_params();


%% Pipeline
% 1) simulate and extract one burst per sample of I
N = numel(S.I.values);
raw = ml.simulate_and_extract_bursts(S.I.values, p0, S.sim);

% 2) arc-length resample curves (Cs) + velocities
[Cs, vels, arc] = curves.resample_bursts_arclength(raw, S.curve.M);

% 3) validate (optional): require equal spike counts
if S.validation.requireSameSpikeCount
    ml.validate_spike_counts(raw);
end

% 4) Geometric Fréchet mean with shift-alignment along arc length
[meanC, aligned, alignedVels, fmInfo] = frechet.mean_shift_aligned(Cs, vels, S.fm);

% 5) reconstruct time parametrizations from speed along the curves
[timeInfo, meanTimeInfo] = curves.reconstruct_time_from_speed(aligned, alignedVels, raw.periods);

% Build mean-curve cumulative time using mean geometry + average speed profile
[meanCumTime, meanPosWrap] = curves.cumtime_from_curve_speed(meanC, meanTimeInfo.meanSp, meanTimeInfo.meanPeriod);
meanTimeInfo.cumTime = meanCumTime;
meanTimeInfo.posWrap = meanPosWrap;

%% Metrics
% Geometric Fréchet mean metrics analysis
metrics = analysis.compute_frechet_metrics(aligned, meanC, arc, fmInfo, S.metrics);
analysis.export_metrics_table(metrics, S.export);

%% Collect outputs
out = struct();
out.settings = S;
out.params = p0;
out.raw = raw;
out.curves = Cs;
out.vels = vels;
out.arc = arc;
out.meanC = meanC;
out.aligned = aligned;
out.alignedVels = alignedVels;
out.metrics = metrics;
out.fmInfo = fmInfo;
out.timeInfo = timeInfo;
out.meanTimeInfo = meanTimeInfo;

%% Plots

% Plot style setting
utils.apply_style(S.plot);
colors = utils.make_colors(N, S.plot.useColor);
 
% Plotting metrics
if S.plot.metrics
    plotting.plot_metrics_summary(metrics, meanC, S.export);
end


% Plotting 2D trajectories
if S.plot.trajectories2D
    plotting.plot_timeseries_2D(raw, meanC, meanTimeInfo, colors, S.export);
    plotting.plot_speed_profiles_2D(timeInfo, meanTimeInfo, raw, colors, S.export);
end

% Plotting 3D trajectories
if S.plot.trajectories3D
    h = plotting.plot_trajectories_3D(aligned, meanC, raw, colors, S.export);
    if S.export.snapshot
        plotting.export_snapshots_3D(h.fig, h.trackers, timeInfo, meanTimeInfo, S.export);
    end
    plotting.animate_trackers_3D(h.trackers, timeInfo, meanTimeInfo, S.anim);
end

end
