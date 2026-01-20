function handles = plot_timeseries(raw, meanC, meanTimeInfo, colors, plotOpts)
%PLOT_TIMESERIES Plot V(t), w(t), y(t) for each burst and overlay mean.
%
% This keeps plotting out of the algorithm. The only "logic" here is mapping
% data to a plot, and optionally exporting.

arguments
    raw (1,1) struct
    meanC (:,3) double
    meanTimeInfo (1,1) struct
    colors (:,3) double
    plotOpts struct = struct()
end

if ~isfield(plotOpts,"export"), plotOpts.export = false; end
if ~isfield(plotOpts,"outDir"),  plotOpts.outDir  = "figures"; end

if plotOpts.export
    if ~isfolder(plotOpts.outDir)
        mkdir(plotOpts.outDir);
    end
end

% Mean curve in time (wrapped)
Cm_w = [meanC; meanC(1,:)];
tmean = meanTimeInfo.cumTime(1:end-1);

names = { 'V', 'w', 'y' };
cols = [1 2 3];
fileTags = { '2D_V', '2D_w', '2D_y' };

titles = { 'V(t) over one burst', 'w(t) over one burst', 'y(t) over one burst' };

ylabels = { 'V', 'w', 'y' };

handles = struct();
handles.fig = gobjects(3,1);

for i = 1:3
    handles.fig(i) = figure('Color','w','Position',[100 100 1400 600]);
    hold on;

    for k = 1:numel(raw.segT)
        tloc = raw.segT{k} - raw.segT{k}(1);
        z = raw.segZ{k}(:, cols(i));
        h = plot(tloc, z, 'Color', colors(k,:), 'DisplayName', sprintf('I=%.4f', raw.I_values(k)));
        % if alpha supported for line objects
        try
            h.Color(4) = 0.7;
        catch
        end
    end

    plot(tmean, Cm_w(1:end-1, cols(i)), 'k', 'LineWidth', 1.5, 'DisplayName', 'mean');
    xlabel('t');
    ylabel(ylabels{i});
    title(titles{i});

    if plotOpts.export
        exportgraphics(handles.fig(i), fullfile(plotOpts.outDir, fileTags{i} + ".pdf"), 'ContentType','vector');
    end
end

end
