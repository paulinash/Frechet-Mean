function metrics = compute_frechet_metrics(curvesAligned, meanC, arc, fmInfo, opts)
%COMPUTE_FRECHET_METRICS Quantitative sanity checks for Fréchet mean results.
%
% Inputs
%   curvesAligned : cell{N} each d x M curve, shift-aligned
%   meanC         : d x M mean curve
%   arc           : struct or numeric info from resampling (optional)
%   fmInfo        : info returned by frechet.mean_shift_aligned (optional)
%   opts          : struct with fields (optional)
%
% Output
%   metrics : struct with distances, variance, modality, curvature/smoothness, etc.

if nargin < 5, opts = struct(); end
opts = fill_defaults(opts);

N = numel(curvesAligned);

% --- Distances to mean ---
d_to_mean = zeros(N,1);
for i = 1:N
    d_to_mean(i) = analysis.curve_distance(curvesAligned{i}, meanC);
end

% Fréchet variance (empirical): mean squared distance to mean
frechetVar = mean(d_to_mean.^2);

% --- Pairwise distances (subsample if large) ---
np = N*(N-1)/2;
pairwiseD = zeros(np,1);
pairs = zeros(np,2);

dopts = struct();
dopts.allowShift = opts.pairwise.allowShift;
dopts.tauGrid    = opts.pairwise.tauGrid;
dopts.interp     = opts.pairwise.interp;

k = 0;
for i = 1:N
    for j = i+1:N
        k = k + 1;
        pairs(k,:) = [i j];
        pairwiseD(k) = analysis.curve_distance(curvesAligned{i}, curvesAligned{j}, dopts);
    end
end

% Compare scale: average pairwise distance vs mean distance to mean
avgPairwise = mean(pairwiseD);
avgToMean   = mean(d_to_mean);

% A common identity in Euclidean settings: mean pairwise squared distance
% relates to variance. Not exact for arbitrary metrics, but still useful
% as a consistency check:
avgPairwiseSq = mean(pairwiseD.^2);

% --- Modality check of d_to_mean distribution ---
modality = analysis.modality_from_hist(d_to_mean, opts.modality);

% --- Curvature / smoothness metrics (mean vs samples) ---
[curvMean, curvSamples] = analysis.curvature_metrics(meanC, curvesAligned, opts.curvature);

% --- Optional: FM iteration diagnostics if available ---
fmDiag = struct();
if nargin >= 4 && ~isempty(fmInfo)
    fmDiag = fmInfo;
end

metrics = struct();
metrics.N = N;

metrics.dist = struct();
metrics.dist.toMean = d_to_mean;
metrics.dist.mean = avgToMean;
metrics.dist.var = var(d_to_mean);
metrics.dist.frechetVar = frechetVar;

metrics.pairwise = struct();
metrics.pairwise.d = pairwiseD;
metrics.pairwise.pairs = pairs;
metrics.pairwise.mean = avgPairwise;
metrics.pairwise.meanSq = avgPairwiseSq;

metrics.sanity = struct();
metrics.sanity.ratio_meanToMeanPairwise = avgToMean / avgPairwise;  % scale check
metrics.sanity.ratio_varToPairwiseSq    = frechetVar / avgPairwiseSq;

metrics.modality = modality;

metrics.curvature = struct();
metrics.curvature.mean = curvMean;
metrics.curvature.samples = curvSamples;

metrics.fm = fmDiag;

% optional: arc info
metrics.arc = arc;

end



function opts = fill_defaults(opts)

% --- Pairwise distance options ---
if ~isfield(opts,'pairwise'), opts.pairwise = struct(); end
if ~isfield(opts.pairwise,'allowShift'), opts.pairwise.allowShift = true; end
if ~isfield(opts.pairwise,'tauGrid'),    opts.pairwise.tauGrid = []; end
if ~isfield(opts.pairwise,'interp'),     opts.pairwise.interp = "pchip"; end

% If tauGrid not provided, we can try to steal it from fmInfo (if present)
% but compute_frechet_metrics doesn't know fmInfo here in fill_defaults.
% So we only default to empty; analysis.curve_distance will then use linspace(0,1,M).

% --- Modality settings ---
if ~isfield(opts,'modality'), opts.modality = struct(); end
if ~isfield(opts.modality,'numBins'),     opts.modality.numBins = 30; end
if ~isfield(opts.modality,'smoothSigma'), opts.modality.smoothSigma = 1.0; end

% --- Curvature settings ---
if ~isfield(opts,'curvature'), opts.curvature = struct(); end
if ~isfield(opts.curvature,'eps'), opts.curvature.eps = 1e-12; end

end
