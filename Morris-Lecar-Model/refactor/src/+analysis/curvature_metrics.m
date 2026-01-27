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
% C: M x d, assumed (approximately) arc-length parameterized.

assert(isnumeric(C) && ismatrix(C), 'C must be a numeric MxD matrix.');
M = size(C,1);
assert(M >= 4, 'Need at least 4 samples to compute second-difference energy of kappa.');

% First differences along the curve parameter (rows)
dC = diff(C, 1, 1);                     % (M-1) x d
ds = sqrt(sum(dC.^2, 2)) + eps0;        % (M-1) x 1

% Unit tangents
T = dC ./ ds;                           % (M-1) x d

% Tangent differences
dT = diff(T, 1, 1);                     % (M-2) x d
ds_mid = (ds(1:end-1) + ds(2:end)) / 2; % (M-2) x 1

% Curvature
kappa = sqrt(sum(dT.^2, 2)) ./ (ds_mid + eps0);  % (M-2) x 1

% Summary stats
m = struct();
m.kappa = kappa;
m.kappa_mean = mean(kappa);
m.kappa_rms  = sqrt(mean(kappa.^2));

% Riemann-sum style integrals
m.total_curvature = sum(kappa .* ds_mid);        % approx ∫ κ ds
m.bending_energy  = sum((kappa.^2) .* ds_mid);   % approx ∫ κ^2 ds

% --- NEW: wobble / roughness of kappa(s) ---
dk = diff(kappa);                   % (M-3) x 1, discrete derivative of kappa vs index
m.kappa_total_variation = sum(abs(dk));

ddk = diff(kappa, 2);               % (M-4) x 1, second difference of kappa
m.kappa_second_difference_energy = sum(ddk.^2);
end


