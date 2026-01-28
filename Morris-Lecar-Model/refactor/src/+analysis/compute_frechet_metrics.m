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

%% Distance: aligned curves to mean curve 
d_to_mean = zeros(N,1);
for i = 1:N
    d_to_mean(i) = analysis.curve_distance(curvesAligned{i}, meanC);
end

%% Fréchet variance (empirical): mean squared distance to mean
frechetVar = mean(d_to_mean.^2);

% opts for Frechet medoid and pairwise distances
dopts = struct();
dopts.allowShift = opts.pairwise.allowShift;
dopts.tauGrid    = opts.pairwise.tauGrid;
dopts.interp     = opts.pairwise.interp;


%% Pairwise distances and Frechet medoid
D = nan(N,N);          %  distance matrix
D(1:N+1:end) = 0;      % diagonal = 0

% nested helper function
function dij = getDist(i,j)
    if i == j
        dij = 0;
        return;
    end
    if i > j
        tmp=i; i=j; j=tmp;
    end
    if ~isnan(D(i,j))
        dij = D(i,j);
        return;
    end
    dij = analysis.curve_distance(curvesAligned{i}, curvesAligned{j}, dopts);
    D(i,j) = dij;
    D(j,i) = dij;
end

% Pairwise distance
fprintf('Computing Pairwise distances of sample curves ');
npTotal = N*(N-1)/2;
maxPairs = 200;


if npTotal <= maxPairs
    % --- compute all pairs ---
    fprintf('using all %d pairs. \n', npTotal);

    pairwiseD = zeros(npTotal,1);
    pairs = zeros(npTotal,2);

    k = 0;
    for i = 1:N
        for j = i+1:N
            k = k + 1;
            pairs(k,:) = [i j];
            pairwiseD(k) = getDist(i,j);
        end
    end
else
    % --- randomly subsample maxPairs unordered pairs ---
    fprintf('sampling %d pairs instead of %d pairs. \n',maxPairs, npTotal);

    pairwiseD = zeros(maxPairs,1);
    pairs = zeros(maxPairs,2);

    k = 0;
    while k < maxPairs
        i = randi(N); j = randi(N);
        if i == j
            continue;
        end
        if i > j
            tmp = i; i = j; j = tmp;
        end
        % avoid duplicates
        if k > 0 && any(pairs(1:k,1) == i & pairs(1:k,2) == j)
            continue;
        end
        k = k + 1;
        pairs(k,:) = [i j];
        pairwiseD(k) = getDist(i,j);
    end
end
avgPairwise = mean(pairwiseD.^2);



% --- Fréchet medoid results (store to locals first) ---
fprintf('Computing Fréchet medoid ');
medoidIdx = NaN;
medoidC = [];
bestCost = NaN;
candIdx = [];
candidateCosts = []; % K random candidates

if npTotal <= maxPairs
    % We computed ALL pairs -> D is complete -> exact medoid
    fprintf('using all %d pairs. \n', npTotal);
    medoidCost = mean(D.^2, 1);           % 1xN
    [bestCost, medoidIdx] = min(medoidCost);
    medoidC = curvesAligned{medoidIdx};
else
    % We subsampled pairs -> D incomplete -> approximate medoid via random candidates
    K = min(N, opts.medoid.numCandidates);   % K samples candidates
    fprintf('by subsampling. \n');

    candIdx = randperm(N, K); % index of sampled candidates

    candidateCosts = zeros(K,1);
    for k = 1:K
        j = candIdx(k);
        s = 0;
        for i = 1:N
            dij = getDist(i,j);      
            s = s + dij^2;
        end
        candidateCosts(k) = s / N;
    end

    [bestCost, kBest] = min(candidateCosts);
    medoidIdx = candIdx(kBest);
    medoidC   = curvesAligned{medoidIdx};
end

d_medoid_to_mean = analysis.curve_distance(medoidC, meanC, dopts);
if frechetVar > 0
    varRatio_medoid = bestCost / frechetVar;
else
    varRatio_medoid = NaN;
end



%% --- Modality check of d_to_mean distribution ---
fprintf('Computing modality stats...\n');
modality = analysis.modality_from_hist(d_to_mean, opts.modality);


%% --- Curvature / smoothness metrics (mean vs samples) ---
% correspondence parametrization ( for plotting curvature only!)
fprintf('Computing curvature stats...\n');
[curvMean_corr, curvSamples_corr] = analysis.curvature_metrics(meanC, curvesAligned, opts.curvature);

% arc length parametrization (for metrics use)
M = size(meanC, 1);                    

meanC_al = resample_uniform_arclength(meanC, M, opts.curvature.eps);
curvesAligned_al = curvesAligned;
for i = 1:N
    curvesAligned_al{i} = resample_uniform_arclength(curvesAligned{i}, M, opts.curvature.eps);
end
[curvMean_al, curvSamples_al] = analysis.curvature_metrics(meanC_al, curvesAligned_al, opts.curvature);


% --- Optional: FM iteration diagnostics if available ---
fmDiag = struct();
if nargin >= 4 && ~isempty(fmInfo)
    fmDiag = fmInfo;
end

metrics = struct();
metrics.N = N;

metrics.dist = struct();
metrics.dist.toMean = d_to_mean;
metrics.dist.frechetVar = frechetVar;

metrics.pairwise = struct();
metrics.pairwise.d = pairwiseD;
metrics.pairwise.mean = avgPairwise;
metrics.pairwise.ratio_meanToMeanPairwise = frechetVar / avgPairwise;  % scale check

metrics.medoid = struct();
metrics.medoid.idx = medoidIdx;         % index of frechet medoid
metrics.medoid.curve = medoidC;         % Frechet medoid
metrics.medoid.cost = bestCost;         % corresponds to frechet Variance for medoid
metrics.medoid.d_to_mean = d_medoid_to_mean; % distance between frechet medoid and frechet mean
metrics.medoid.varRatio = varRatio_medoid; % ratio between FV(medoid) and FV(mean)

if ~isempty(candIdx)
    metrics.medoid.candidates = candIdx;
    metrics.medoid.candidateCosts = candidateCosts;
end

metrics.modality = modality;

metrics.curvature = struct();
metrics.curvature.corr   = struct('mean', curvMean_corr, 'samples', curvSamples_corr);
metrics.curvature.arclen = struct('mean', curvMean_al,   'samples', curvSamples_al);

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

% --- Medoid settings ---
if ~isfield(opts,'medoid'), opts.medoid = struct(); end
if ~isfield(opts.medoid,'numCandidates'), opts.medoid.numCandidates = 20; end


% --- Curvature settings ---
if ~isfield(opts,'curvature'), opts.curvature = struct(); end
if ~isfield(opts.curvature,'eps'), opts.curvature.eps = 1e-12; end

end



function C_u = resample_uniform_arclength(C, Mnew, eps0)
%RESAMPLE_UNIFORM_ARCLENGTH_Md Resample an M x d curve to uniform arc length.

assert(isnumeric(C) && ismatrix(C), 'C must be numeric MxD.');
M = size(C,1);
assert(M >= 2, 'Need at least 2 points to resample.');

dC = diff(C, 1, 1);                      % (M-1) x d
ds = sqrt(sum(dC.^2, 2)) + eps0;         % (M-1) x 1
s  = [0; cumsum(ds)];                    % M x 1
L  = s(end);

if L < eps0
    C_u = repmat(C(1,:), Mnew, 1);       % Mnew x d
    return;
end

s_new = linspace(0, L, Mnew).';
C_u = interp1(s, C, s_new, 'pchip');     % Mnew x d
end


