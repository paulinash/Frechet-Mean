function export_snapshots_3D(fig3D, trackers, timeInfo, meanTimeInfo, plotOpts)
%EXPORT_SNAPSHOTS_3D Export 3D snapshots with trackers at selected times.
%
% Inputs
%   fig3D        : figure handle
%   trackers     : struct with fields .sample (Nx1) and .mean
%   timeInfo     : struct with fields .cumTime{k} (M+1x1), .posWrap{k} (M+1x3)
%   meanTimeInfo : struct with fields .cumTime (M+1x1), .posWrap (M+1x3), .meanPeriod
%   plotOpts.outDir (char/string) output directory
%   plotOpts.snapshotN (int) number of snapshots (default 9)
%
% Notes
%   - Snapshot times are chosen uniformly over one mean period, excluding endpoint.
%   - Each sample tracker is advanced by its own reconstructed period (mod).

arguments
    fig3D (1,1) matlab.ui.Figure
    trackers (1,1) struct
    timeInfo (1,1) struct
    meanTimeInfo (1,1) struct
    plotOpts struct = struct()
end

if ~isfield(plotOpts,'outDir'), plotOpts.outDir = "figures"; end
if ~isfield(plotOpts,'snapshotN'), plotOpts.snapshotN = 9; end
plotOpts.outDir = string(plotOpts.outDir);

outDir = plotOpts.outDir;
if isstring(outDir), outDir = char(outDir); end
if ~isfolder(outDir)
    mkdir(outDir);
end

nFrames = plotOpts.snapshotN;
meanPeriod = meanTimeInfo.meanPeriod;

frameTimes = linspace(0, meanPeriod, nFrames+1);
frameTimes(end) = [];

N = numel(timeInfo.cumTime);

for f = 1:numel(frameTimes)
    t_global = frameTimes(f);

    % sample trackers
    for k = 1:N
        cumt = timeInfo.cumTime{k};
        posW = timeInfo.posWrap{k};
        tloc = mod(t_global, cumt(end));
        px = interp1(cumt, posW(:,1), tloc, 'pchip');
        py = interp1(cumt, posW(:,2), tloc, 'pchip');
        pz = interp1(cumt, posW(:,3), tloc, 'pchip');
        set(trackers.sample(k), 'XData', px, 'YData', py, 'ZData', pz);
    end

    % mean tracker
    cumt = meanTimeInfo.cumTime;
    posW = meanTimeInfo.posWrap;
    tloc = mod(t_global, cumt(end));
    px = interp1(cumt, posW(:,1), tloc, 'pchip');
    py = interp1(cumt, posW(:,2), tloc, 'pchip');
    pz = interp1(cumt, posW(:,3), tloc, 'pchip');
    set(trackers.mean, 'XData', px, 'YData', py, 'ZData', pz);

    drawnow;

    exportgraphics(fig3D, fullfile(outDir, sprintf('3D_frame_%02d.pdf', f)), ...
        'Resolution', 600, 'ContentType', 'image');
end

end