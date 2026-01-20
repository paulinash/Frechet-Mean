function h = plot_trajectories_3d(curves, meanC, raw, colors, plotOpts)
%PLOT_TRAJECTORIES_3D 3D trajectories and trackers.
%
% Returns a struct with fields:
%   h.fig
%   h.trackers.sample (Nx1)
%   h.trackers.mean

arguments
    curves (1,:) cell
    meanC (:,3) double
    raw (1,1) struct
    colors (:,3) double
    plotOpts struct = struct()
end

if ~isfield(plotOpts,"export"), plotOpts.export = false; end
if ~isfield(plotOpts,"outDir"),  plotOpts.outDir  = "figures"; end
plotOpts.outDir = string(plotOpts.outDir);

if plotOpts.export && ~isfolder(plotOpts.outDir)
    mkdir(plotOpts.outDir);
end

N = numel(curves);

fig3D = figure('Color','w','Position',[200 100 900 700]);
hold on; grid on;

for k = 1:N
    C = curves{k};
    hline = plot3([C(:,1); C(1,1)], [C(:,2); C(1,2)], [C(:,3); C(1,3)], ...
        'Color', colors(k,:), 'LineWidth', 1.2, 'DisplayName', sprintf('I=%.4f', raw.I_values(k)));
    try
        hline.Color(4) = 0.3;
    catch
    end
end

plot3([meanC(:,1);meanC(1,1)], [meanC(:,2);meanC(1,2)], [meanC(:,3);meanC(1,3)], ...
    'k', 'LineWidth', 3, 'DisplayName', 'mean curve');

xlabel('V'); ylabel('w'); zlabel('y');
if isfield(raw,'spikeCounts') && ~isempty(raw.spikeCounts)
    title(sprintf('3D trajectories with %d spikes', raw.spikeCounts(1)));
else
    title('3D trajectories');
end
view([-30 20]);

if plotOpts.export
    exportgraphics(fig3D, fullfile(plotOpts.outDir, "3D.pdf"), 'ContentType','image', 'Resolution', 600);
end

% Trackers
trackers = struct();
trackers.sample = gobjects(N,1);
for k=1:N
    C = curves{k};
    trackers.sample(k) = plot3(C(1,1), C(1,2), C(1,3), 'o', ...
        'Color', colors(k,:), 'MarkerFaceColor', colors(k,:), 'MarkerSize', 7, 'HandleVisibility','off');
end
trackers.mean = plot3(meanC(1,1), meanC(1,2), meanC(1,3), 'ok', ...
    'MarkerSize', 10, 'MarkerFaceColor', 'k', 'DisplayName', 'mean tracker');

h = struct();
h.fig = fig3D;
h.trackers = trackers;
end
