function vdp_frechet_mean_tracking_closestpoint_tracker()
% Compute Fréchet mean using closest-point method and animate trackers

clear; close all; clc;

%% Parameters
epsilon = 0.1;                  
a_vals = linspace(-0.5, 0.9, 2); 
Tfinal = 400;                   
dt = 0.01;                      
M = 500;                        % number of points to resample each curve
useArcLength = true;            

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
    
    % Uniform time sampling
    if useArcLength
        curves{k} = resampleArcLength(z_ss, M);
    else
        curves{k} = interp1(linspace(0,1,size(z_ss,1)), z_ss, linspace(0,1,M));
    end
end

%% Compute "closest-point" Fréchet mean
curve1 = curves{1};
curve2 = curves{2};
meanCurve = zeros(size(curve1));
closestIdxCurve2 = zeros(size(curve1,1),1); % store indices for tracker

for i = 1:size(curve1,1)
    p1 = curve1(i,:);
    dists = sqrt(sum((curve2 - p1).^2,2));
    [~, idxMin] = min(dists);
    p2 = curve2(idxMin,:);
    closestIdxCurve2(i) = idxMin;
    
    meanCurve(i,:) = (p1 + p2)/2;
end

%% Animation to track curves step by step
figure; hold on; grid on;
colors = lines(length(a_vals));
axis equal
xlabel('x'); ylabel('y');
title('Van der Pol solutions with closest-point Fréchet mean (tracked points)');

% Plot full trajectories
plot(curve1(:,1), curve1(:,2), '--','Color',colors(1,:));
plot(curve2(:,1), curve2(:,2), '--','Color',colors(2,:));
plot(meanCurve(:,1), meanCurve(:,2), 'k--');

% Initialize moving points
points = gobjects(3,1);
points(1) = plot(curve1(1,1), curve1(1,2), 'o', 'MarkerFaceColor', colors(1,:), ...
                 'MarkerEdgeColor','k', 'MarkerSize',8);
points(2) = plot(curve2(closestIdxCurve2(1),1), curve2(closestIdxCurve2(1),2), 'o', ...
                 'MarkerFaceColor', colors(2,:), 'MarkerEdgeColor','k', 'MarkerSize',8);
points(3) = plot(meanCurve(1,1), meanCurve(1,2), 'o', 'MarkerFaceColor', 'k', ...
                 'MarkerEdgeColor','k','MarkerSize',8);

% Animate
for i = 1:size(curve1,1)
    % Curve 1 moves normally
    set(points(1),'XData',curve1(i,1),'YData',curve1(i,2));
    % Curve 2 jumps to closest point used for mean
    idx2 = closestIdxCurve2(i);
    set(points(2),'XData',curve2(idx2,1),'YData',curve2(idx2,2));
    % Mean point
    set(points(3),'XData',meanCurve(i,1),'YData',meanCurve(i,2));
    
    drawnow;
    pause(0.1);
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