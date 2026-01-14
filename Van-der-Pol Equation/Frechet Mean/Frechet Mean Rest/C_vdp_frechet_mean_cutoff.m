function vdp_frechet_mean_tracking_cutoff()
% Visualize two Van der Pol solutions with Fréchet mean and adjustable cutoff
% The first curve (blue) can start at a later index, e.g., index 50.

clear; close all; clc;

%% Parameters
epsilon = 0.1;
a_vals = linspace(-0.9, 0.9, 2);
Tfinal = 400;
dt = 0.01;
M = 500;
useArcLength = false;

%% -------------------- CUT-OFF FOR BLUE CURVE --------------------
cutoffIdx = 75;  % manually choose starting index for first curve (blue)

%% Preallocate curves
curves = cell(length(a_vals),1);

%% Solve Van der Pol for each 'a'
for k = 1:length(a_vals)
    a = a_vals(k);
    vdp = @(t,z) [ z(2) - z(1)^3/3 + z(1); epsilon*(a - z(1)) ];
    z0 = [-1;0];
    [t,z] = ode45(vdp,0:dt:Tfinal,z0);
    
    % Use last half of trajectory
    z_ss = z(t >= Tfinal/2,:);
    
    % Resample along arc length or uniform time
    if useArcLength
        curves{k} = resampleArcLength(z_ss, M);
    else
        curves{k} = interp1(linspace(0,1,size(z_ss,1)), z_ss, linspace(0,1,M));
    end
end

%% Apply cutoff to first curve
curves{1} = curves{1}(cutoffIdx:end, :);

% Shorten second curve to match first for mean calculation
lenMin = size(curves{1},1);
curves{2} = curves{2}(1:lenMin, :);

%% Compute Fréchet mean (pointwise average)
allCurves = zeros(lenMin,2,length(curves));
for k = 1:length(curves)
    allCurves(:,:,k) = curves{k};
end
meanCurve = mean(allCurves,3);

%% ------------------- Animation -------------------
figure; hold on; grid on;
colors = lines(length(a_vals));
axis equal
xlabel('x'); ylabel('y');
title('Van der Pol solutions with cutoff and Fréchet mean');

% Plot entire trajectories in faint colors
for k = 1:length(a_vals)
    plot(curves{k}(:,1), curves{k}(:,2), '--','Color',colors(k,:));
end
plot(meanCurve(:,1), meanCurve(:,2),'k--');

% Initialize moving points
points = gobjects(length(a_vals)+1,1);
for k = 1:length(a_vals)
    points(k) = plot(curves{k}(1,1), curves{k}(1,2), 'o', ...
        'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor','k','MarkerSize',8);
end
points(end) = plot(meanCurve(1,1), meanCurve(1,2), 'o', ...
    'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',8);

% Animate step by step
for i = 1:lenMin
    for k = 1:length(a_vals)
        set(points(k),'XData',curves{k}(i,1),'YData',curves{k}(i,2));
    end
    set(points(end),'XData',meanCurve(i,1),'YData',meanCurve(i,2));
    drawnow;
    pause(1);
end

%% --------- Helper function: Arc-length resampling ---------
function curveOut = resampleArcLength(curveIn, M)
    diffs = diff(curveIn,1,1);
    segLen = sqrt(sum(diffs.^2,2));
    cumLen = [0; cumsum(segLen)];
    t_new = linspace(0, cumLen(end), M);
    x_new = interp1(cumLen, curveIn(:,1), t_new);
    y_new = interp1(cumLen, curveIn(:,2), t_new);
    curveOut = [x_new(:), y_new(:)];
end

end
