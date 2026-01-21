function [meanC, aligned, alignedVels, info] = mean_shift_aligned(curves, vels, opts)
%MEAN_SHIFT_ALIGNED Shift-align periodic curves and compute their Fréchet mean.
%
%   [meanC, aligned, alignedVels, info] = frechet.mean_shift_aligned(curves, vels, opts)
%
% Inputs
%   curves{k}      : MxD array, uniform samples in normalized arc length s \in [0,1)
%   vels{k}        : MxD array, velocity vectors (w.r.t. original time) resampled to same s-grid
%   opts.maxIter   : maximum iterations
%   opts.tol       : relative convergence tolerance on the mean curve
%   opts.tauGrid   : vector of candidate shifts in [0,1]
%
% Outputs
%   meanC          : MxD Fréchet mean (pointwise Euclidean mean after best shifts)
%   aligned{k}     : MxD aligned curves
%   alignedVels{k} : MxD aligned velocities
%   info           : diagnostics (taus, errors, history)
%
% Notes
% - This is a *circular shift along the arc-length parameter* (a rotation of the
%   parametrization), not a full reparameterization.

arguments
    curves (1,:) cell
    vels   (1,:) cell
    opts struct = struct()
end

% ---- defaults (robust library style) ----
if ~isfield(opts,"maxIter"), opts.maxIter = 40; end
if ~isfield(opts,"tol"),     opts.tol = 1e-8; end
if ~isfield(opts,"tauGrid"), opts.tauGrid = []; end
if ~isfield(opts,"verbose"), opts.verbose = false; end
% ----------------------------------------

N = numel(curves);
M = size(curves{1},1);
D = size(curves{1},2);

% Basic validation
for k = 1:N
    if size(curves{k},1) ~= M || size(curves{k},2) ~= D
        error('All curves must have the same size MxD.');
    end
    if size(vels{k},1) ~= M || size(vels{k},2) ~= D
        error('All velocity arrays must match curve sizes MxD.');
    end
end

if isempty(opts.tauGrid)
    % Default: coarse grid is risky; use M points like your original
    opts.tauGrid = linspace(0,1,M).';
else
    opts.tauGrid = opts.tauGrid(:);
end

% Arc grid (implicit): assume curves are sampled at s_i = i/M, i=0..M-1.
sqFull = linspace(0,1,M+1).';
sq = sqFull(1:end-1);

% Initialization: pointwise average
meanC = mean(cat(3, curves{:}), 3);
aligned = curves;
alignedVels = vels;

info = struct();
info.taus = zeros(N, opts.maxIter);
info.errors = zeros(N, opts.maxIter);
info.meanChange = zeros(opts.maxIter,1);

for iter = 1:opts.maxIter % Iterate until convergence
    mean_old = meanC;
    
    % Iterate over each curve
    for k = 1:N
        bestErr = Inf; bestTau = 0; bestCurve = [];
        
        % Iterate over each possible shift
        for tau = opts.tauGrid.'

            % Rotate each curve 
            sq_shifted = mod(sq + tau, 1);

            % Interpolate curve onto shifted grid (periodic; extrap ok because mod)
            gammaShift = interp1(sq, curves{k}, sq_shifted, 'pchip', 'extrap');

            % Pointwise distance to currenct mean
            err = mean(sum((gammaShift - meanC).^2, 2));
            if err < bestErr
                bestErr = err;
                bestTau = tau;
                bestCurve = gammaShift;
            end
        end

        % Choose best alignment to current mean    
        aligned{k} = bestCurve;
        sq_shifted_best = mod(sq + bestTau, 1);
        alignedVels{k} = interp1(sq, vels{k}, sq_shifted_best, 'pchip', 'extrap');

        info.taus(k, iter) = bestTau;
        info.errors(k, iter) = bestErr;

        fprintf('iter %d/%d, curve %d: tau=%.4f, err=%.6g\n', iter, opts.maxIter, k, bestTau, bestErr);
    end

    % Update mean using mean of optimal aligned curves
    meanC = mean(cat(3, aligned{:}), 3);

    % Relative change of old to new mean
    relChange = norm(meanC(:) - mean_old(:)) / (norm(mean_old(:)) + eps);
    info.meanChange(iter) = relChange;

    if relChange < opts.tol
        info.nIter = iter;
        info.taus = info.taus(:,1:iter);
        info.errors = info.errors(:,1:iter);
        info.meanChange = info.meanChange(1:iter);
        fprintf('Fréchet mean converged in %d iterations (relChange=%.3g).\n', iter, relChange);
        return;
    end
end

info.nIter = opts.maxIter;
fprintf('Fréchet mean reached maxIter=%d (last relChange=%.3g).\n', opts.maxIter, info.meanChange(end));

end
