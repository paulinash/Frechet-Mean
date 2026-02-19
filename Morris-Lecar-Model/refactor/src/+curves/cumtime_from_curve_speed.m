function [cumTime, posWrap] = cumtime_from_curve_speed(curve, speedVec, period)
%CUMTIME_FROM_CURVE_SPEED Build a time parametrization from a speed profile.
%
% We assume the curve is sampled uniformly in normalized arc length and is
% periodic. Using local increments:
%   dt_i = ds_i / speed_i
% where ds_i is the Euclidean distance between consecutive samples on the wrapped
% curve. Then we rescale so that cumTime(end) = period.
%
% Inputs
%   curve    : MxD (not wrapped)
%   speedVec : Mx1 (speed magnitudes at each sample)
%   period   : scalar
%
% Outputs
%   cumTime  : (M+1)x1 cumulative time
%   posWrap  : (M+1)xD curve with the first point appended at the end

M = size(curve,1);
if numel(speedVec) ~= M
    error('speedVec must have length M (number of curve samples).');
end

% Wrap curve
posWrap = [curve; curve(1,:)];

% Arc length increments
dC = diff(posWrap, 1, 1);
ds = sqrt(sum(dC.^2, 2));

% Speed magnitudes
sp = speedVec(:);
sp(sp < 1e-12) = 1e-12;

% Local time increments
dt_local = ds ./ sp;
cumTime = [0; cumsum(dt_local)];

% Rescale to match requested period
cumTime = cumTime * (period / cumTime(end));

end
