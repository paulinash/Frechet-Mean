function S = config_exp_ml_frechet_mean(projectRoot)
%CONFIG_EXP_ML_FRECHET_MEAN Settings for exp_ml_frechet_mean.
%
% Usage:
%   projectRoot = fileparts(fileparts(mfilename('fullpath')));
%   S = config_exp_ml_frechet_mean(projectRoot);

arguments
    projectRoot (1,1) string
end

S = struct();

%% I sweep
S.I.low = 44.8;
S.I.high = 44.9;
S.I.N = 5;
S.I.values = linspace(S.I.low, S.I.high, S.I.N);

%% Simulation
S.sim.Tfinal = 8000;
S.sim.dt = 0.05;
S.sim.tspan = 0:S.sim.dt:S.sim.Tfinal;
S.sim.z0 = [-20; 0; 6];
S.sim.transientFraction = 0.5; % drop first 50% to obtain non transient

%% Curve representation
S.curve.M = 2000; % points along arc length

%% Fréchet mean
S.fm.maxIter = 20;
S.fm.tol = 1e-8;
S.fm.tauGrid = linspace(0,1,S.curve.M).';  % column vector is often nicer

%% Validation
S.validation.requireSameSpikeCount = false;

%% Plot styling / global figure behavior
S.plot.useColor = true;

S.plot.style = "paper";     % "paper" or "screen"
S.plot.fontSize = 12;
S.plot.lineWidth = 1.25;
S.plot.figureColor = "w";
S.plot.useLatex = true;
S.plot.metrics = true;
S.plot.trajectories2D = false;
S.plot.trajectories3D = true;

%% Export settings
S.export.export = true;
S.export.snapshot = false;
S.export.outDir = fullfile(projectRoot,'figures');
S.export.snapshotN = 9;

%% Animation
S.anim.dt_play = 0.02;
S.anim.speed_factor = 20;
S.anim.maxSeconds = 40;

%% Metrics
S.metrics = struct();
S.metrics.pairwise = struct();
S.metrics.pairwise.allowShift = true;
S.metrics.pairwise.tauGrid = S.fm.tauGrid;     % reuse the same grid
S.metrics.pairwise.interp = "pchip";
S.metrics.modality.numBins = 30;
S.metrics.modality.smoothSigma = 0.0;

end
