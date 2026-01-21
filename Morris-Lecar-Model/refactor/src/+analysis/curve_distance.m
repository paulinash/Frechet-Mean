function [d, info] = curve_distance(C1, C2, opts)
%CURVE_DISTANCE Distance consistent with frechet.mean_shift_aligned objective.
%
%   d = analysis.curve_distance(C1, C2)
%   d = analysis.curve_distance(C1, C2, opts)
%   [d, info] = analysis.curve_distance(...)
%
% C1, C2: MxD arrays sampled uniformly in normalized arc-length s in [0,1).
%
% Default distance (no shift search):
%   d = sqrt( mean( ||C1(s_i) - C2(s_i)||^2 ) )
% which is exactly sqrt of the per-curve error used in mean_shift_aligned.
%
% If opts.allowShift = true, then performs shift-invariant distance:
%   d = min_tau sqrt( mean( ||C1(s_i+tau) - C2(s_i)||^2 ) )
% using the same periodic pchip interpolation as mean_shift_aligned.
%
% opts fields (all optional):
%   allowShift : false by default
%   tauGrid    : vector of candidate shifts in [0,1]. If empty and allowShift,
%                uses linspace(0,1,M).
%   interp     : 'pchip' (default), should match mean_shift_aligned
%   verbose    : false
%
% info returns:
%   info.bestTau, info.bestErr, info.errs (if shift search)

arguments
    C1 (:,:) double
    C2 (:,:) double
    opts struct = struct()
end

% defaults
if ~isfield(opts,'allowShift'), opts.allowShift = false; end
if ~isfield(opts,'tauGrid'),    opts.tauGrid = []; end
if ~isfield(opts,'interp'),     opts.interp = 'pchip'; end
if ~isfield(opts,'verbose'),    opts.verbose = false; end

assert(all(size(C1) == size(C2)), 'C1 and C2 must be same size MxD.');
M = size(C1,1);

info = struct();

% ---- base (no shift) ----
if ~opts.allowShift
    % match objective: mean(sum((C1-C2).^2,2)) then sqrt
    err = mean(sum((C1 - C2).^2, 2));
    d = sqrt(err);
    info.bestTau = 0;
    info.bestErr = err;
    return;
end

% ---- shift-invariant (min over tauGrid) ----
if isempty(opts.tauGrid)
    tauGrid = linspace(0,1,M).';
else
    tauGrid = opts.tauGrid(:);
end

% implicit arc grid like in mean_shift_aligned
sqFull = linspace(0,1,M+1).';
sq = sqFull(1:end-1);

bestErr = Inf;
bestTau = 0;

errs = zeros(numel(tauGrid),1);

for t = 1:numel(tauGrid)
    tau = tauGrid(t);

    sq_shifted = mod(sq + tau, 1);
    C1shift = interp1(sq, C1, sq_shifted, opts.interp, 'extrap');

    err = mean(sum((C1shift - C2).^2, 2));
    errs(t) = err;

    if err < bestErr
        bestErr = err;
        bestTau = tau;
    end
end

d = sqrt(bestErr);

info.bestTau = bestTau;
info.bestErr = bestErr;
info.errs = errs;

if opts.verbose
    fprintf('curve_distance: bestTau=%.4f, bestErr=%.6g, d=%.6g\n', bestTau, bestErr, d);
end

end
