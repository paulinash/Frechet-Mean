function animate_trackers_3d(trackers, timeInfo, meanTimeInfo, plotOpts)
%ANIMATE_TRACKERS_3D Live animation of trackers moving along curves.
%
% This is intentionally isolated because it is a "long-running" loop and often
% the first thing you want to comment out when doing analysis or exporting.
%
% Usage:
%   plotting.animate_trackers_3d(trackers, timeInfo, meanTimeInfo, plotOpts)
%
% plotOpts fields (all optional):
%   .dt_play (default 0.02)
%   .speed_factor (default 20)
%
% Stop with Ctrl+C.

if nargin < 4
    plotOpts = struct();
end
if ~isfield(plotOpts,'dt_play'), plotOpts.dt_play = 0.02; end
if ~isfield(plotOpts,'speed_factor'), plotOpts.speed_factor = 20; end

fprintf('Animation running... Press Ctrl+C to stop.\n');

N = numel(timeInfo.cumTime);
t_global = 0;

dt_play = plotOpts.dt_play;
speed_factor = plotOpts.speed_factor;

while true
    t_global = t_global + speed_factor*dt_play;

    for k = 1:N
        cumt = timeInfo.cumTime{k};
        posW = timeInfo.posWrap{k};
        tloc = mod(t_global, cumt(end));
        px = interp1(cumt, posW(:,1), tloc, 'pchip');
        py = interp1(cumt, posW(:,2), tloc, 'pchip');
        pz = interp1(cumt, posW(:,3), tloc, 'pchip');
        set(trackers.sample(k), 'XData', px, 'YData', py, 'ZData', pz);
    end

    cumt = meanTimeInfo.cumTime;
    posW = meanTimeInfo.posWrap;
    tloc = mod(t_global, cumt(end));
    px = interp1(cumt, posW(:,1), tloc, 'pchip');
    py = interp1(cumt, posW(:,2), tloc, 'pchip');
    pz = interp1(cumt, posW(:,3), tloc, 'pchip');
    set(trackers.mean, 'XData', px, 'YData', py, 'ZData', pz);

    drawnow limitrate;
    pause(dt_play);
end

end
