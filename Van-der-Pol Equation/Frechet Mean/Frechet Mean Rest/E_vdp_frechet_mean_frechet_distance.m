function vdp_frechet_mean_starting_at_intersection()
% Two Van der Pol limit cycles
% Align them by choosing a common intersection as the starting point
% Then compute Fréchet matching + mean curve with no jumping.

clear; close all; clc;

%% Parameters
epsilon = 0.1;
a_vals = [-0.5, 0.9];
Tfinal = 400;
dt = 0.01;
M = 500;
useArcLength = false;

%% Generate the curves
curves = cell(2,1);

for k = 1:2
    a = a_vals(k);
    vdp = @(t,z)[ z(2) - z(1)^3/3 + z(1) ; epsilon*(a - z(1)) ];
    [t,z] = ode45(vdp, 0:dt:Tfinal, [-1;0]);
    z_ss = z(t >= Tfinal/2, :);

    if useArcLength
        curves{k} = resampleArcLength(z_ss, M);
    else
        curves{k} = interp1(linspace(0,1,size(z_ss,1)), z_ss, linspace(0,1,M));
    end
end

curve1 = curves{1};
curve2 = curves{2};

%% --------------------------------------------------------------------
%   1. FIND INTERSECTION AND ALIGN CURVES BY THAT START POINT
%% --------------------------------------------------------------------
fprintf("Finding intersection point...\n");

[idx1_inter, idx2_inter] = findCurveIntersection(curve1, curve2);

% Rotate curves so that the intersection is at index 1
curve1 = circshift(curve1, -(idx1_inter-1), 1);
curve2 = circshift(curve2, -(idx2_inter-1), 1);

fprintf("Intersection found at curve1(%d) and curve2(%d).\n", ...
    idx1_inter, idx2_inter);
fprintf("Curves rotated to start at intersection.\n");

%% --------------------------------------------------------------------
%   2. FRECHET MATCHING
%% --------------------------------------------------------------------
fprintf("Computing discrete Fréchet alignment...\n");

path = discreteFrechetPath(curve1, curve2);

idx1 = path(:,1);
idx2 = path(:,2);

% Resample the path back to M points
ii = linspace(1, length(idx1), M);
idx1i = round(interp1(1:length(idx1), idx1, ii));
idx2i = round(interp1(1:length(idx2), idx2, ii));

meanCurve = 0.5*(curve1(idx1i,:) + curve2(idx2i,:));

fprintf("Done.\n");

%% --------------------------------------------------------------------
%   3. ANIMATION
%% --------------------------------------------------------------------

figure; hold on; grid on; axis equal;
colors = lines(2);
xlabel('x'); ylabel('y');
title('Fréchet Matching Starting at Intersection');

plot(curve1(:,1), curve1(:,2), '--', 'Color', colors(1,:));
plot(curve2(:,1), curve2(:,2), '--', 'Color', colors(2,:));
plot(meanCurve(:,1), meanCurve(:,2), 'k--');

p1 = plot(curve1(idx1i(1),1), curve1(idx1i(1),2), 'o', ...
    'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor','k', 'MarkerSize',8);
p2 = plot(curve2(idx2i(1),1), curve2(idx2i(1),2), 'o', ...
    'MarkerFaceColor', colors(2,:), 'MarkerEdgeColor','k', 'MarkerSize',8);
pm = plot(meanCurve(1,1), meanCurve(1,2), 'o', ...
    'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',8);

for k = 1:M
    set(p1,'XData',curve1(idx1i(k),1),'YData',curve1(idx1i(k),2));
    set(p2,'XData',curve2(idx2i(k),1),'YData',curve2(idx2i(k),2));
    set(pm,'XData',meanCurve(k,1),'YData',meanCurve(k,2));
    drawnow;
    pause(0.05);
end


%% ====================================================================
%                       Helper Functions
%% ====================================================================

function curveOut = resampleArcLength(curveIn, M)
    diffs = diff(curveIn,1,1);
    segLen = sqrt(sum(diffs.^2,2));
    cumLen = [0; cumsum(segLen)];
    t_new = linspace(0, cumLen(end), M);
    curveOut = [ interp1(cumLen, curveIn(:,1), t_new).' ...
                 interp1(cumLen, curveIn(:,2), t_new).' ];
end

%% ---- Find approximate intersection of two curves (closest pair) ----
function [i1, i2] = findCurveIntersection(C1, C2)
    D = pdist2(C1, C2);         % M × M distance matrix
    [i1, i2] = find(D == min(D(:)), 1);
end

%% ---- Discrete Fréchet Matching ----
function path = discreteFrechetPath(A, B)

    M = size(A,1); N = size(B,1);
    D = sqrt( (A(:,1)-B(:,1)').^2 + (A(:,2)-B(:,2)').^2 );
    F = zeros(M,N);

    F(1,1) = D(1,1);
    for i = 2:M, F(i,1) = max(F(i-1,1),D(i,1)); end
    for j = 2:N, F(1,j) = max(F(1,j-1),D(1,j)); end

    for i = 2:M
        for j = 2:N
            F(i,j) = max(D(i,j), min([F(i-1,j), F(i-1,j-1), F(i,j-1)]));
        end
    end

    i=M; j=N;
    path = [i,j];
    while ~(i==1 && j==1)
        if i==1
            j=j-1;
        elseif j==1
            i=i-1;
        else
            [~,k] = min([F(i-1,j), F(i-1,j-1), F(i,j-1)]);
            if k==1, i=i-1;
            elseif k==2, i=i-1; j=j-1;
            else, j=j-1;
            end
        end
        path = [i,j; path];
    end
end

end