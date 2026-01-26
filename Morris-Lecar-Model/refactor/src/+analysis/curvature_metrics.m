function [meanMetrics, sampleMetrics] = curvature_metrics(meanC, curves, opts)
%CURVATURE_METRICS Curvature/smoothness comparison for mean and samples.
%
% Inputs
%   meanC  : M x d mean curve (samples along rows)
%   curves : cell{N} of M x d sample curves
%   opts.eps : small epsilon for numerical stability
%
% Outputs
%   meanMetrics   : struct with curvature stats for meanC
%   sampleMetrics : struct array (N x 1) with curvature stats for each curve

eps0 = opts.eps;

meanMetrics = curvature_of_curve(meanC, eps0);

N = numel(curves);

% Preallocate struct array with correct fields
sampleMetrics = repmat(curvature_of_curve(curves{1}, eps0), N, 1);

% Fill all entries (including i=1 for clarity/robustness)
for i = 1:N
    sampleMetrics(i) = curvature_of_curve(curves{i}, eps0);
end

end


function m = curvature_of_curve(C, eps0)
% C: M x d, assumed (approximately) arc-length parametrized in the sense
% that successive point spacing is roughly uniform along the curve.

assert(isnumeric(C) && ismatrix(C), 'C must be a numeric MxD matrix.');
M = size(C,1);
d = size(C,2);
assert(M >= 3, 'Need at least 3 samples along the curve to compute curvature.');

% First differences along the curve parameter (rows)
dC = diff(C, 1, 1);                     % (M-1) x d
ds = sqrt(sum(dC.^2, 2)) + eps0;        % (M-1) x 1, segment lengths

% Unit tangent vectors for each segment
T = dC ./ ds;                           % (M-1) x d

% Differences of tangents
dT = diff(T, 1, 1);                     % (M-2) x d
ds_mid = (ds(1:end-1) + ds(2:end)) / 2; % (M-2) x 1

% Discrete curvature: ||dT|| / ds
kappa = sqrt(sum(dT.^2, 2)) ./ (ds_mid + eps0);  % (M-2) x 1

% Summary stats
m = struct();
m.kappa = kappa;                         % (M-2) x 1
m.kappa_mean = mean(kappa);              % average curvature
m.kappa_rms  = sqrt(mean(kappa.^2));     % RMS curvature

% Discrete integrals (Riemann sums)
m.total_curvature = sum(kappa .* ds_mid);        % approx ∫ κ ds
m.bending_energy  = sum((kappa.^2) .* ds_mid);   % approx ∫ κ^2 ds

% Second difference smoothness proxy (parameterization-dependent)
ddC = diff(C, 2, 1);                     % (M-2) x d
m.second_diff_rms = sqrt(mean(sum(ddC.^2, 2)));

end
