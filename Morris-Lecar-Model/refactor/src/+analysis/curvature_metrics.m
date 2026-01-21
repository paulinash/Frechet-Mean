function [meanMetrics, sampleMetrics] = curvature_metrics(meanC, curves, opts)
%CURVATURE_METRICS Curvature/smoothness comparison for mean and samples.
eps0 = opts.eps;

meanMetrics = curvature_of_curve(meanC, eps0);

N = numel(curves);
% build once to get the correct field layout, then preallocate consistently
sampleMetrics = curvature_of_curve(curves{1}, eps0);
sampleMetrics = repmat(sampleMetrics, N, 1);

for i = 2:N
    sampleMetrics(i) = curvature_of_curve(curves{i}, eps0);
end

end



function m = curvature_of_curve(C, eps0)
% C: d x M, assumed arc-length parameter (approximately)
dC = diff(C,1,2);                 % d x (M-1)
ds = sqrt(sum(dC.^2,1)) + eps0;   % 1 x (M-1)

T = dC ./ ds;                     % unit tangent approx, d x (M-1)

dT = diff(T,1,2);                 % d x (M-2)
ds_mid = (ds(1:end-1) + ds(2:end))/2;  % 1 x (M-2)

kappa = sqrt(sum(dT.^2,1)) ./ (ds_mid + eps0);   % 1 x (M-2)

% Summary stats
m = struct();
m.kappa = kappa;
m.kappa_mean = mean(kappa);
m.kappa_rms  = sqrt(mean(kappa.^2));
m.kappa_max  = max(kappa);

% Discrete integrals (approx)
m.total_curvature = sum(kappa .* ds_mid);
m.bending_energy  = sum((kappa.^2) .* ds_mid);

% Another "smoothness" proxy: second difference norm (geometry-only)
ddC = diff(C,2,2);                           % d x (M-2)
m.second_diff_rms = sqrt(mean(sum(ddC.^2,1)));

end
