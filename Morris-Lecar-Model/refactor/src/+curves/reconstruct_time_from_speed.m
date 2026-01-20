function [timeInfo, meanTimeInfo] = reconstruct_time_from_speed(curves, vels, periods)
%RECONSTRUCT_TIME_FROM_SPEED Reconstruct a time parameterization from speed.
%
% For each curve sampled in arc length (periodic), we use
%   dt_local = ds / ||v||
% along the discrete arc grid. We then rescale to match the detected burst period.
%
% Inputs:
%   curves{k}: Mx3 (not wrapped)
%   vels{k}:   Mx3 (speed vectors along same arc grid)
%   periods:   Nx1 burst durations
%
% Outputs:
%   timeInfo.cumTime{k}: (M+1)x1 cumulative time, wrapped (includes endpoint)
%   timeInfo.posWrap{k}: (M+1)x3 curve wrapped with first point appended
%   timeInfo.speeds(k,:): (1xM) speed magnitudes
%   meanTimeInfo has analogous fields for the average speed profile and average period

N = numel(curves);
M = size(curves{1},1);

timeInfo = struct();
timeInfo.cumTime = cell(N,1);
timeInfo.posWrap = cell(N,1);
timeInfo.speeds = zeros(N,M);

for k = 1:N
    C = curves{k};
    Vv = vels{k};

    Cw = [C; C(1,:)];
    dC = diff(Cw,1,1);
    ds = sqrt(sum(dC.^2,2));

    sp = sqrt(sum(Vv.^2,2));
    sp(sp < 1e-12) = 1e-12;
    timeInfo.speeds(k,:) = sp.';

    dt_local = ds ./ sp;
    cumt = [0; cumsum(dt_local)];
    cumt = cumt * (periods(k) / cumt(end));

    timeInfo.cumTime{k} = cumt;
    timeInfo.posWrap{k} = Cw;
end

% Mean time parametrization: average speed profile + mean period
meanSp = mean(timeInfo.speeds,1).';
meanSp(meanSp < 1e-12) = 1e-12;
meanPeriod = mean(periods);

meanTimeInfo = struct();
meanTimeInfo.meanSp = meanSp;
meanTimeInfo.meanPeriod = meanPeriod;
meanTimeInfo.cumTime = [];
meanTimeInfo.posWrap = [];

% We don't have mean curve here; caller can pass it to plotting. The mean cumTime
% depends on the mean curve geometry, so compute it in plotting using mean curve.
% However, for convenience, we return mean speed and mean period.
end
