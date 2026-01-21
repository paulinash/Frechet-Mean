function out = exp_ml_frechet_mean()
%EXP_ML_FRECHET_MEAN Morris–Lecar bursting + (shift-aligned) Fréchet mean.
%
% This is an experiment runner: keep it short. Put algorithms in src/.
%
clear; close all; clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(projectRoot,'src')));

%% Settings (equivalent to your original script, but grouped)
S = struct();
S.I.low = 44.8;
S.I.high = 45.2;
S.I.N = 10;
S.I.values = linspace(S.I.low, S.I.high, S.I.N);

S.sim.Tfinal = 8000;
S.sim.dt = 0.05;
S.sim.tspan = 0:S.sim.dt:S.sim.Tfinal;
S.sim.z0 = [-20; 0; 6];
S.sim.transientFraction = 0.5; % drop first 50% to obtain non transient

S.curve.M = 2000; % points along arc length

S.fm.maxIter = 20;
S.fm.tol = 1e-8;
S.fm.tauGrid = linspace(0,1,S.curve.M);

S.validation.requireSameSpikeCount = false;

% Plot styling / global figure behavior
S.plot.useColor = true;

% Global plot style (paper)
S.plot.style = "paper";     % "paper" or "screen"
S.plot.fontSize = 12;       % good default for papers
S.plot.lineWidth = 1.25;
S.plot.figureColor = "w";
S.plot.useLatex = true;     % set interpreters to latex

% Export settings (used by plotting/export functions)
S.export.export = true;
S.export.snapshot = false;
S.export.outDir = fullfile(projectRoot,'figures');
S.export.snapshotN = 9;           % for export_snapshots_3d (optional)

% Animation settings (used only by animate_trackers_3d)
S.anim.dt_play = 0.02;
S.anim.speed_factor = 20;
S.anim.maxSeconds = 40;

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

% 4) Fréchet mean with shift-alignment along arc length
[meanC, aligned, alignedVels, fmInfo] = frechet.mean_shift_aligned(Cs, vels, S.fm);

% 5) reconstruct time parametrizations from speed along the curves
[timeInfo, meanTimeInfo] = curves.reconstruct_time_from_speed(aligned, alignedVels, raw.periods);

% Build mean-curve cumulative time using mean geometry + average speed profile
[meanCumTime, meanPosWrap] = curves.cumtime_from_curve_speed(meanC, meanTimeInfo.meanSp, meanTimeInfo.meanPeriod);
meanTimeInfo.cumTime = meanCumTime;
meanTimeInfo.posWrap = meanPosWrap;


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
out.fmInfo = fmInfo;
out.timeInfo = timeInfo;
out.meanTimeInfo = meanTimeInfo;

%% Plots (all plotting separated)
utils.apply_style(S.plot);
colors = utils.make_colors(N, S.plot.useColor);


plotting.plot_timeseries_2D(raw, meanC, meanTimeInfo, colors, S.export);
plotting.plot_speed_profiles_2D(timeInfo, meanTimeInfo, raw, colors, S.export);

h = plotting.plot_trajectories_3D(aligned, meanC, raw, colors, S.export);

if S.export.snapshot
    plotting.export_snapshots_3D(h.fig, h.trackers, timeInfo, meanTimeInfo, S.export);
end

plotting.animate_trackers_3D(h.trackers, timeInfo, meanTimeInfo, S.anim);
end
