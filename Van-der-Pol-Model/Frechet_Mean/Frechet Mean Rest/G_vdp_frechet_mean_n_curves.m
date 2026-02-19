function vdp_frechet_mean_with_tracking
% Iterative Fréchet (Karcher) mean for N Van-der-Pol periodic curves
% With moving tracking points (animated markers).

clear; close all; clc;

%% USER OPTIONS
paramMode = 'arc';     % 'arc' or 'time'
epsilon   = 0.1;

a_vals = linspace(-0.8,0.95,7);    % ANY number of curves
Na = numel(a_vals);

Tfinal = 400;
dt = 0.005;
M = 300;

maxIter = 30;
tol     = 1e-6;

%% ================================================================
%   SIMULATE + EXTRACT ONE PERIOD + RESAMPLE
%% ================================================================
curves_res = cell(Na,1);

for k = 1:Na
    a = a_vals(k);
    vdp = @(t,z) [z(2) - z(1)^3/3 + z(1) ;
                   epsilon*(a - z(1)) ];

    tspan = 0:dt:Tfinal;
    z0 = [-1; 0];

    [t,z] = ode45(vdp, tspan, z0);
    
    % only consider non transient phase
    idx = t >= Tfinal/2;
    z2 = z(idx,:);
    t2 = t(idx);

    % find peaks and period of solution
    x = z2(:,1);
    [~,locs] = findpeaks(x,'MinPeakProminence',0.05,'MinPeakDistance',10);
    if numel(locs) < 2
        error("Not enough peaks for a = %.3f", a);
    end
    i1 = locs(end-1);
    i2 = locs(end);

    seg  = z2(i1:i2,:);
    tseg = t2(i1:i2);

    % chose arclength or time parametrization with M points
    switch paramMode
        case 'arc'
            curves_res{k} = resampleArcLengthClosed(seg, M);
        case 'time'
            tU = linspace(tseg(1), tseg(end), M);
            curves_res{k} = interp1(tseg, seg, tU);
    end
end

%% ================================================================
%   ITERATIVE FRÉCHET / KARCHER MEAN
%% ================================================================
aligned = curves_res;

sGrid = linspace(0,1,M+1);
sGrid = sGrid(1:end-1);   % periodic grid
tauGrid = linspace(0,1,M);

% Initial mean = average of unaligned curves
meanCurve = zeros(M,2);
for k = 1:Na
    meanCurve = meanCurve + curves_res{k};
end
meanCurve = meanCurve / Na;

for iter = 1:maxIter
    mean_old = meanCurve;

    % --- Align each curve to current mean ---
    for k = 1:Na
        C = curves_res{k};

        bestScore = Inf;
        bestTau   = 0;

        for tau = tauGrid
            sShift = mod(sGrid + tau, 1);

            % continuous reparametrization
            Csh = interp1(sGrid, C, sShift, 'pchip');

            diff2 = sum((Csh - meanCurve).^2, 2);
            score = mean(diff2);

            if score < bestScore
                bestScore = score;
                bestTau   = tau;
            end
        end

        % apply optimal continuous shift
        sShift = mod(sGrid + bestTau, 1);
        aligned{k} = interp1(sGrid, C, sShift, 'pchip');

    end

    % --- Update mean ---
    meanCurve = zeros(M,2);
    for k = 1:Na
        meanCurve = meanCurve + aligned{k};
    end
    meanCurve = meanCurve / Na;

    % --- Check convergence ---
    relChange = norm(meanCurve(:) - mean_old(:)) / (norm(mean_old(:))+eps);
    fprintf("Iter %d: rel change = %.3e\n",iter,relChange);

    if relChange < tol
        fprintf("Converged after %d iterations.\n", iter);
        break
    end
end

%% ================================================================
%   PLOT CURVES + TRACKING POINTS ANIMATION
%% ================================================================
figure('Color','w','Position',[100 100 1000 700]);
hold on; axis equal; grid on;
if strcmp(paramMode,'arc')
    title('Fréchet Mean of Van der Pol Cycles using arc-length parametrization');
else
    title('Fréchet Mean of Van der Pol Cycles using time parametrization');
end


colors = lines(Na);

% ---- Plot aligned curves ----
for k = 1:Na
    plot(aligned{k}(:,1), aligned{k}(:,2), 'Color',colors(k,:), ...
        'LineWidth',1.2, 'DisplayName',sprintf('a=%.2f',a_vals(k)));
end

% ---- Plot mean curve ----
plot(meanCurve(:,1), meanCurve(:,2),'k','LineWidth',2,'DisplayName','mean curve');

% ---- Plot centers ----
for k = 1:Na
    C = aligned{k};
    cen = mean(C,1);
    plot(cen(1), cen(2), 's', 'Color',colors(k,:), ...
        'MarkerFaceColor',colors(k,:), 'MarkerSize',9, 'HandleVisibility','off');
end
cenMean = mean(meanCurve,1);
plot(cenMean(1), cenMean(2), 'ks','MarkerFaceColor','k','MarkerSize',10, 'DisplayName','center');


legend('Location','southeast');
xlabel('x'); ylabel('y');


%% ================================================================
%   TRACKING POINTS (MOVING MARKERS)
%% ================================================================
% Create a marker for each curve
hCurve = gobjects(Na,1);
for k = 1:Na
    C = aligned{k};
    hCurve(k) = plot(C(1,1), C(1,2),'o', ...
        'Color',colors(k,:), 'MarkerFaceColor',colors(k,:), ...
        'MarkerSize',8, 'HandleVisibility','off');
end

% Marker for mean curve
hMean = plot(meanCurve(1,1), meanCurve(1,2), 'ko', ...
             'MarkerFaceColor','k', 'MarkerSize',9, 'DisplayName','tracker');

%% Animation loop
while true
    for i = 1:M
        for k = 1:Na
            Ci = aligned{k}(i,:);
            set(hCurve(k),'XData',Ci(1),'YData',Ci(2));
        end

        set(hMean,'XData',meanCurve(i,1),'YData',meanCurve(i,2));

        drawnow limitrate;
        pause(0.02);
    end
end

end

%% ================================================================
% Helper: arc-length parametrization
%% ================================================================
function C = resampleArcLengthClosed(seg, M)
    seg2 = [seg; seg(1,:)];

    d = diff(seg2,1,1);
    segLen = sqrt(sum(d.^2,2));
    cumLen = [0; cumsum(segLen)];
    L = cumLen(end);

    s = linspace(0,L,M+1);
    x = interp1(cumLen, seg2(:,1), s);
    y = interp1(cumLen, seg2(:,2), s);

    C = [x(1:end-1).'  y(1:end-1).'];
end