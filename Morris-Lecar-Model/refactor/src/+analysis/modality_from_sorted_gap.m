function out = modality_from_sorted_gap(x, opts)
%MODALITY_FROM_SORTED_GAP Detects "pillars" via biggest jump in sorted distances.
%
% Idea:
%   Sort distances: xs = sort(x).
%   Compute jumps:  dx = diff(xs).
%   Look for large *relative* jump: dx / (xs + eps).
%   A big jump means two (or more) groups: close-to-mean vs far-from-mean.
%
% opts fields (simple defaults):
%   opts.gapThreshold     : relative jump threshold to declare split (default 1.0)
%   opts.minGroupSize     : minimum points per group (default 2)

x = x(:);
x = x(isfinite(x));
n = numel(x);

out = struct();
out.n = n;

if n < 2
    out.class = "single_group";
    out.sorted = x;
    out.splitIndex = [];
    out.gapScore = 0;
    out.relGap = [];
    return;
end

xs = sort(x);
dx = diff(xs);

% relative gap: robust across scales (0.05..0.5 vs 0.5..12)
eps0 = 1e-12;
relGap = dx ./ (xs(1:end-1) + eps0);

% biggest relative jump
[gapScore, k] = max(relGap);

% defaults
if nargin < 2, opts = struct(); end
if ~isfield(opts,'gapThreshold'), opts.gapThreshold = 1.0; end
if ~isfield(opts,'minGroupSize'), opts.minGroupSize = 2; end

% avoid splits that isolate 1 point
isValidSplit = (k >= opts.minGroupSize) && ((n-k) >= opts.minGroupSize);

isTwoGroup = isValidSplit && (gapScore >= opts.gapThreshold);

out.sorted = xs;
out.diff = dx;
out.relGap = relGap;
out.gapScore = gapScore;
out.splitIndex = k;
out.splitValueLeft  = xs(k);
out.splitValueRight = xs(k+1);

if isTwoGroup
    out.class = "two_or_more_groups";
else
    out.class = "single_group";
end
end
