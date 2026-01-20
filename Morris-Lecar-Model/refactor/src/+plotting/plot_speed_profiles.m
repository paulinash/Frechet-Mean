function hfig = plot_speed_profiles(timeInfo, meanTimeInfo, raw, colors, plotOpts)
%PLOT_SPEED_PROFILES Plot reconstructed speed profiles vs reconstructed time.

arguments
    timeInfo (1,1) struct
    meanTimeInfo (1,1) struct
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

hfig = figure('Color','w','Position',[100 100 1400 600]);
hold on;

N = size(timeInfo.speeds,1);
for k = 1:N
    tk = timeInfo.cumTime{k}(1:end-1);
    h = plot(tk, timeInfo.speeds(k,:), 'Color', colors(k,:), 'DisplayName', sprintf('I=%.4f', raw.I_values(k)));
    try
        h.Color(4) = 0.35;
    catch
    end
end

% Mean speed (already computed in curves.reconstruct_time_from_speed)
plot(meanTimeInfo.cumTime(1:end-1), meanTimeInfo.meanSp, 'k', 'LineWidth', 2.0, 'DisplayName', 'mean speed');

xlabel('reconstructed time t (one burst)');
ylabel('speed ||z\_dot||');
title('Speed profiles vs reconstructed time');

if plotOpts.export
    exportgraphics(hfig, fullfile(plotOpts.outDir, "speed_vs_time.pdf"), 'ContentType','vector');
end
end
