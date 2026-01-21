function apply_style(plotOpts)
%APPLY_STYLE Set global plotting defaults (fonts, interpreters, linewidths).
%
% Call this once at the start of an experiment, before creating figures.
%
% plotOpts fields (optional):
%   .style        "paper" (default) or "screen"
%   .fontSize     default axes font size
%   .lineWidth    default line width
%   .figureColor  figure background color
%   .useLatex     true/false set latex interpreters

if nargin < 1 || isempty(plotOpts), plotOpts = struct(); end

if ~isfield(plotOpts,'style'),       plotOpts.style = "paper"; end
if ~isfield(plotOpts,'fontSize'),    plotOpts.fontSize = 12; end
if ~isfield(plotOpts,'lineWidth'),   plotOpts.lineWidth = 1.25; end
if ~isfield(plotOpts,'figureColor'), plotOpts.figureColor = "w"; end
if ~isfield(plotOpts,'useLatex'),    plotOpts.useLatex = true; end

% ---- global defaults for this MATLAB session (or until changed) ----
set(groot, ...
    'defaultFigureColor', plotOpts.figureColor, ...
    'defaultAxesFontSize', plotOpts.fontSize, ...
    'defaultAxesLineWidth', 1.0, ...
    'defaultLineLineWidth', plotOpts.lineWidth);

if plotOpts.useLatex
    set(groot, ...
        'defaultTextInterpreter', 'latex', ...
        'defaultAxesTickLabelInterpreter', 'latex', ...
        'defaultLegendInterpreter', 'latex');
else
    set(groot, ...
        'defaultTextInterpreter', 'tex', ...
        'defaultAxesTickLabelInterpreter', 'tex', ...
        'defaultLegendInterpreter', 'tex');
end

% Optional tweaks for different modes
switch string(plotOpts.style)
    case "paper"
        % reasonable paper defaults already set above
    case "screen"
        % example: bigger for presentations
        set(groot, 'defaultAxesFontSize', max(plotOpts.fontSize, 14));
    otherwise
        % leave as-is
end
end
