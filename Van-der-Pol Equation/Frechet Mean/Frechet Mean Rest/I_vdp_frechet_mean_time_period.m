function vdp_frechet_mean_with_tracking
% Iterative Fréchet (Karcher) mean for N Van-der-Pol periodic curves
% With moving tracking points (animated markers)
% CHANGED: now computes mean period when paramMode = 'time'

clear; close all; clc;

%% USER OPTIONS
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
periods    = zeros(Na,1);   % ADDED: store periods for 'time' mode

for k = 1:Na
    a = a_vals(k);
    vdp = @(t,z) [ z(2) - z(1)^3/3 + z(1) ;
                   epsilon*(a - z(1)) ];

    tspan = 0:dt:Tfinal;
    z0 = [-1; 0];

    [t,z] = ode45(vdp, tspan, z0);

    idx = t >= Tfinal/2;
    z2 = z(idx,:);
    t2 = t(idx);

    x = z2(:,1);
    [~,locs] = findpeaks(x,'MinPeakProminence',0.05,'MinPeakDistance',10);
    if numel(locs) < 2
        error("Not enough peaks for a = %.3f", a);
    end

    i1 = locs(end-1);
    i2 = locs(end);

    seg  = z2(i1:i2,:);
    tseg = t2(i1:i2);

    % ADDED: store period for this sample
    periods(k) = tseg(end) - tseg(1);
    fprintf('Period of curve %d: %2.4f\n', k, periods(k));

    
    % CHANGED: normalize time to [0,1] for interpolation but keep period info
    tU = linspace(tseg(1), tseg(end), M);
    curves_res{k} = interp1(tseg, seg, tU);
   
end

%% ================================================================
%   ITERATIVE FRÉCHET / KARCHER MEAN
%% ================================================================
aligned = curves_res;
sGrid  = linspace(0,1,M+1);
sGrid  = sGrid(1:end-1);   % periodic normalized time grid
tauGrid = linspace(0,1,M);   % phase shift search grid

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

            % continuous time reparametrization
            Csh = interp1(sGrid, C, sShift, 'pchip');

            diff2 = sum((Csh - meanCurve).^2, 2);
            score = mean(diff2);

            if score < bestScore
                bestScore = score;
                bestTau   = tau;
            end
        end

        % apply optimal continuous phase shift
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
%   MEAN PERIOD 
%% ================================================================

meanPeriod = mean(periods);   % ADDED: average the periods
fprintf('Mean period: %2.4f\n', meanPeriod);


%% ================================================================
%   PLOT CURVES + TRACKING POINTS ANIMATION
%% ================================================================
figure('Color','w','Position',[100 100 1000 700]);
hold on; axis equal; grid on;
title('Fréchet Mean of Van der Pol Cycles using time parametrization');


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
hCurve = gobjects(Na,1);
for k = 1:Na
    C = aligned{k};
    hCurve(k) = plot(C(1,1), C(1,2),'o', ...
        'Color',colors(k,:), 'MarkerFaceColor',colors(k,:), ...
        'MarkerSize',8, 'HandleVisibility','off');
end

hMean = plot(meanCurve(1,1), meanCurve(1,2), 'ko', ...
             'MarkerFaceColor','k', 'MarkerSize',9, 'DisplayName','tracker');

%% Animation loop
% CHANGED: include mean period in time-parametrization
t_global = 0;
dt_play = 0.02;
speed_factor=5;

while true
    t_global = t_global + speed_factor*dt_play;

    for k = 1:Na
        % Determine local phase along curve
        
        Tk = periods(k);    % sample's true period
        
        tloc = mod(t_global, Tk);

        % linear interpolation along curve
        % scale index to [1,M] for resampling
        idxf = (tloc / Tk) * (M-1) + 1;
        i1 = floor(idxf); i2 = i1+1;
        alpha = idxf - i1;

        % handle wrapping
        if i2 > M
            i2 = 1; 
        end

        Ci = (1-alpha)*aligned{k}(i1,:) + alpha*aligned{k}(i2,:);
        set(hCurve(k),'XData',Ci(1),'YData',Ci(2));
    end

    % mean curve animation (scale by mean period)
    Tk_m = meanPeriod;
   
    tloc_m = mod(t_global, Tk_m);
    idxf_m = (tloc_m / Tk_m) * (M-1) + 1;
    i1 = floor(idxf_m); i2 = i1+1;
    alpha = idxf_m - i1;
    if i2 > M, i2 = 1; end
    Cm = (1-alpha)*meanCurve(i1,:) + alpha*meanCurve(i2,:);
    set(hMean,'XData',Cm(1),'YData',Cm(2));

    drawnow limitrate;
    pause(dt_play);
end

end


