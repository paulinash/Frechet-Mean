function vdp_optimal_shift_mean_curve()
% Compute one-period solutions for Van der Pol with parameters a=-0.3 and 0.3
% Optimize the rotational alignment by searching the shift that minimizes
% the global distance to the resulting mean curve.
%
% Includes:
%   - arc-length or time parametrization
%   - global optimal rotational alignment
%   - mean curve construction
%   - center alignment
%   - tracking points animation in both plots

clear; close all; clc;

%% USER OPTION: parametrization mode
paramMode = 'time';   % 'arc' or 'time'

%% Parameters
epsilon = 0.1;
a_vals = [-0.9, 0.9];
Tfinal = 400;
dt = 0.005;
M = 300;     % number of samples per period

%% ================================================================
%   SOLVE VAN DER POL AND EXTRACT ONE PERIOD FOR BOTH PARAMETERS
%% ================================================================
curves_raw = cell(2,1);
tcell = cell(2,1);

for k = 1:2
    a = a_vals(k);
    vdp = @(t,z) [ z(2) - z(1).^3/3 + z(1); epsilon*(a - z(1)) ];
    z0 = [-1; 0];
    tspan = 0:dt:Tfinal;

    [t,z] = ode45(vdp, tspan, z0);
    idx = t >= Tfinal/2;       % use second half only
    curves_raw{k} = z(idx,:);
    tcell{k} = t(idx);
end

%% ================================================================
%   EXTRACT ONE PERIOD USING PEAKS OF x
%% ================================================================
curves_period = cell(2,1);
t_period = cell(2,1);

for k = 1:2
    x = curves_raw{k}(:,1);
    t = tcell{k};

    [~,locs] = findpeaks(x,'MinPeakProminence',0.05,'MinPeakDistance',10);
    if numel(locs) < 2
        error("Not enough peaks for one full period.");
    end

    i1 = locs(end-1); 
    i2 = locs(end);

    seg = curves_raw{k}(i1:i2,:);
    tseg = t(i1:i2);

    curves_period{k} = seg;
    t_period{k} = tseg;
end

%% ================================================================
%   RESAMPLE CURVES ACCORDING TO PARAMETRIZATION MODE
%% ================================================================
curves_resampled = cell(2,1);

switch paramMode
    case 'arc'
        for k = 1:2
            curves_resampled{k} = resampleArcLengthClosed(curves_period{k}, M);
        end
    case 'time'
        for k = 1:2
            tseg = t_period{k};
            t_uniform = linspace(tseg(1), tseg(end), M);
            curves_resampled{k} = interp1(tseg, curves_period{k}, t_uniform);
        end
    otherwise
        error("Invalid mode. Use 'arc' or 'time'.");
end

curve1 = curves_resampled{1};
curve2 = curves_resampled{2};

%% ================================================================
%   OPTIMAL SHIFT: CHOOSE STARTING POINT THAT MINIMIZES GLOBAL MEAN DISTANCE
%% ================================================================
bestShift = 0;
bestScore = Inf;

for kShift = 0:M-1
    curve2_shifted = circshift(curve2, -kShift, 1);
    meanCurve_tmp = 0.5*(curve1 + curve2_shifted);

    dist1 = sqrt(sum((curve1 - meanCurve_tmp).^2, 2));
    dist2 = sqrt(sum((curve2_shifted - meanCurve_tmp).^2, 2));
    score = mean(dist1 + dist2);     % global distance measure

    if score < bestScore
        bestScore = score;
        bestShift = kShift;
    end
end

fprintf("Optimal shift found: %d  (score = %.6f)\n", bestShift, bestScore);

%% Apply best shift
curve2_best = circshift(curve2, -bestShift, 1);
curve1_best = curve1;

%% Final mean curve
meanCurve = 0.5*(curve1_best + curve2_best);


%% ================================================================
%   FIGURE 1: ORIGINAL POSITION, OPTIMAL SHIFT, MEAN CURVE + TRACKING
%% ================================================================
figure('Color','w','Position',[100 100 900 650]); hold on; axis equal; grid on;
title(sprintf('Optimal alignment by global distance (%s parametrization)', paramMode));

plot(curve1_best(:,1), curve1_best(:,2),'LineWidth',1.2,'DisplayName','curve 1');
plot(curve2_best(:,1), curve2_best(:,2),'LineWidth',1.2,'DisplayName','curve 2');
plot(meanCurve(:,1), meanCurve(:,2),'k','LineWidth',2,'DisplayName','mean curve');

h1 = plot(curve1_best(1,1), curve1_best(1,2),'bo','MarkerFaceColor','b');
h2 = plot(curve2_best(1,1), curve2_best(1,2),'ro','MarkerFaceColor','r');
h3 = plot(meanCurve(1,1),meanCurve(1,2),'ko','MarkerFaceColor','k');


% ------------------------------------------------------------------
% NEW: Plot centers of the three curves (original, before center alignment)
% ------------------------------------------------------------------
center1_fig1 = mean(curve1_best,1);
center2_fig1 = mean(curve2_best,1);
centerMean_fig1 = mean(meanCurve,1);

plot(center1_fig1(1), center1_fig1(2), 'bs', 'MarkerFaceColor','b', ...
     'MarkerSize',10, 'DisplayName','center curve 1');

plot(center2_fig1(1), center2_fig1(2), 'rs', 'MarkerFaceColor','r', ...
     'MarkerSize',10, 'DisplayName','center curve 2');

plot(centerMean_fig1(1), centerMean_fig1(2), 'ks', 'MarkerFaceColor','k', ...
     'MarkerSize',10, 'DisplayName','center mean curve');
legend('Location','best');
xlabel('x'); ylabel('y');

%% Tracking animation
for k = 1:M
    set(h1,'XData',curve1_best(k,1),'YData',curve1_best(k,2));
    set(h2,'XData',curve2_best(k,1),'YData',curve2_best(k,2));
    set(h3,'XData',meanCurve(k,1),'YData',meanCurve(k,2));
    drawnow;
end


fprintf("Finished.\n");

%% ================================================================
% Helper function: arc-length parametrization for closed curves
%% ================================================================
function C = resampleArcLengthClosed(seg, M)
    seg2 = [seg; seg(1,:)];

    diffs = diff(seg2,1,1);
    segLen = sqrt(sum(diffs.^2,2));
    cumLen = [0; cumsum(segLen)];
    L = cumLen(end);

    snew = linspace(0,L,M+1);
    x = interp1(cumLen, seg2(:,1), snew);
    y = interp1(cumLen, seg2(:,2), snew);

    C = [x(1:end-1).' y(1:end-1).'];
end

end