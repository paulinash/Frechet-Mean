function vdp_frechet_mean_arclength_dynamic
% Iterative Fréchet (Karcher) mean for N Van-der-Pol periodic curves
% Also extracts periods + velocities and visualizes speed via
% moving markers (trajectories) on each curve and the mean curve.
% Uses an arclength parametrization to preserve geometrical information
% Calculates mean period and speed to preserve dynamical information

clear; close all; clc;

%% USER OPTIONS
epsilon   = 0.1;
dynamics = false;
useColor = true;

% Samples of a
a_vals = linspace(-0.6,0.9,40);
Na = numel(a_vals);

% Simulation time and timesteps for ODE solver
Tfinal = 400;
dt = 0.005;
M = 300;             % geometric samples per curve

% Frechet mean iteration settings
maxIter = 40;
tol     = 1e-8;

%% STORAGE
curves_res = cell(Na,1);
periods = zeros(Na,1);
vel_time = cell(Na,1);
raw_time_samples = cell(Na,1);

%% ================================================================
%   SIMULATE + EXTRACT PERIOD + RESAMPLE + VELOCITIES
%% ================================================================
for k = 1:Na
    % Define Van der Pol System
    a = a_vals(k);
    vdp = @(t,z) [ z(2) - z(1).^3/3 + z(1) ; epsilon*(a - z(1)) ];
    
    % Simulate ODE
    tspan = 0:dt:Tfinal;
    z0 = [-1; 0];       % Starting point of trajectories
    [t,z] = ode45(vdp, tspan, z0);

    % Extract second half of simulation for periodic behavior
    idx = t >= Tfinal/2;
    z2 = z(idx,:);
    t2 = t(idx);
    x = z2(:,1);

    % Find peaks of x(t) to find oscillation periods
    [~,locs] = findpeaks(x,'MinPeakProminence',0.05,'MinPeakDistance',10);
    if numel(locs) < 2
        error("Not enough peaks for a = %.3f", a);
    end
    i1 = locs(end-1); i2 = locs(end);
    
    % Store periodic segments
    seg  = z2(i1:i2,:);
    tseg = t2(i1:i2);

    raw_time_samples{k}.pos = seg;
    raw_time_samples{k}.t   = tseg;

    periods(k) = tseg(end) - tseg(1);

    % get velocities of sample curves
    vx = gradient(seg(:,1), tseg);
    vy = gradient(seg(:,2), tseg);

    % compute cumulative arc-length
    d = diff(seg); segLen = sqrt(sum(d.^2,2));
    cumLen = [0; cumsum(segLen)];
    s_norm = cumLen / cumLen(end);     % normalize arclength to [0,1]

    sq = linspace(0,1,M+1);            % common domain

    x_s = interp1(s_norm, seg(:,1), sq);
    y_s = interp1(s_norm, seg(:,2), sq);
    C = [x_s(1:end-1).' y_s(1:end-1).'];

    vx_s = interp1(s_norm, vx, sq);
    vy_s = interp1(s_norm, vy, sq);
    V_s = [vx_s(1:end-1).' vy_s(1:end-1).'];

    curves_res{k} = C;
    vel_time{k} = V_s;

    % Print curve characteristics
    fprintf('I=%.2f : period dur=%.3f \n', a_vals(k), periods(k));
end
fprintf("Mean period = %.6f\n", mean(periods));

%% ================================================================
%   ITERATIVE FRÉCHET / KARCHER MEAN
%% ================================================================
aligned = curves_res;
alignedVel = vel_time;
tauGrid = linspace(0,1,M);

% Initial guess for mean curve: pointwise average
meanCurve = zeros(M,2);
for k = 1:Na
        meanCurve = meanCurve + curves_res{k}; 
end
meanCurve = meanCurve / Na;

% Iterative Karcher mean
for iter = 1:maxIter
    mean_old = meanCurve;
    
    % Iterate over each curve
    for k = 1:Na
        C = curves_res{k}; bestScore = Inf;

        % Find best alignment to mean curve for each curve
        for tau = tauGrid
            sq_shift = mod(sq(1:end-1) + tau, 1);
            Csh = interp1(sq(1:end-1), C, sq_shift, 'pchip');
            score = mean(sum((Csh - meanCurve).^2,2));
            if score < bestScore
                bestScore = score;
                bestCurve = Csh;
                bestTau = tau;
            end
        end
        aligned{k} = bestCurve;

        % align velocity consistently
        sq_shift = mod(sq(1:end-1) + bestTau, 1);
        alignedVel{k} = interp1(sq(1:end-1), vel_time{k}, sq_shift, 'pchip');
        fprintf('Curve %d best rotation tau = %.4f (error = %.6f)\n', k, bestTau, bestScore);
    end
    
    % Compute new mean
    meanCurve = zeros(M,2);
    for k = 1:Na
        meanCurve = meanCurve + aligned{k}; 
    end
    meanCurve = meanCurve / Na;

    % Check for convergence
    relChange = norm(meanCurve(:)-mean_old(:)) / (norm(mean_old(:))+eps);
    fprintf("Iter %d: rel change = %.3e\n", iter, relChange);
    if relChange < tol
        fprintf("Converged.\n"); 
        break; 
    end
end

%% ================================================================
%   MEAN PERIOD AND MEAN VELOCITY
%% ================================================================
meanPeriod = mean(periods);

% compute speed profiles
allSpeeds = zeros(Na,M);
for k = 1:Na
    V = alignedVel{k};
    allSpeeds(k,:) = sqrt(sum(V.^2,2));
end
meanSpeed = mean(allSpeeds,1).';   % M×1

%% ================================================================
%   PLOT GEOMETRY + MOVING VELOCITY MARKERS
%% ================================================================
% Initialize plot
figure('Color','w','Position',[100 100 1100 700]); hold on; axis equal; grid on;
if dynamics
    title('Fréchet Mean of periodic Van Der Pol solutions with dynamics');
else
    title('Fréchet Mean of periodic Van Der Pol solutions without dynamics');
end

if useColor
    colors = lines(Na);            % one color per sample
    colors = 0.5*colors + 0.5;   % mix with white

else
    gray = [0.75 0.75 0.75];       % light gray
    colors = repmat(gray, Na, 1);  % same gray for all samples
end
% Plot sample and mean curves
for k = 1:Na
    plot([aligned{k}(:,1);  aligned{k}(1,1)], [aligned{k}(:,2);  aligned{k}(1,2)],'Color',colors(k,:),'LineWidth',1.0, 'DisplayName',sprintf('a=%.2f',a_vals(k)));
end
plot([meanCurve(:,1); meanCurve(1,1)], [meanCurve(:,2); meanCurve(1,2)],'k','LineWidth',2,'DisplayName','mean curve');
exportgraphics(gcf,'Figures_vdp/arc_param/fig_frechet_mean.pdf','ContentType','vector');


% Create moving markers showing fast-slow dynamic for sample and mean curves
hCurve = gobjects(Na,1);
for k = 1:Na
    C = aligned{k};
    hCurve(k) = plot(C(1,1), C(1,2),'o','Color',colors(k,:), ...
        'MarkerFaceColor',colors(k,:),'MarkerSize',7,'HandleVisibility','off');
end
hMean = plot(meanCurve(1,1), meanCurve(1,2),'ko','MarkerFaceColor','k','MarkerSize',9,'DisplayName','tracker');

% Label and legend
xlabel('x'); ylabel('y');

%% ================================================================
%   PRECOMPUTE MOTION PARAMETRIZATION (dynamics vs. uniform)
%   Then run ONE animation loop (shared).
%% ================================================================

% Common animation settings
fprintf("Animation running... Press Ctrl+C to stop.\n");
dt_play      = 0.01;
t_global     = 0;

% We'll build, for each curve k:
%   paramGrid{k} : (M+1)x1 increasing parameter (time-like)
%   posGrid{k}   : (M+1)x2 positions along the closed curve
% and for the mean:
%   paramMean : (M+1)x1
%   posMean   : (M+1)x2
paramGrid = cell(Na,1);
posGrid   = cell(Na,1);

% Always wrap the mean curve once
posMean = [meanCurve; meanCurve(1,:)];

if dynamics
    speed_factor = 5;
    % --- Dynamics mode: reconstruct true time from ds / speed, then scale to period ---
    for k = 1:Na
        C = aligned{k};      % Mx2
        V = alignedVel{k};   % Mx2

        % periodic wrap of positions
        Cw = [C; C(1,:)];    % (M+1)x2

        % arc-length increments (M x 1)
        dC = diff(Cw,1,1);
        ds = sqrt(sum(dC.^2,2));

        % speed magnitude (M x 1)
        sp = sqrt(sum(V.^2,2));
        sp(sp < 1e-12) = 1e-12;

        % local time increments (M x 1)
        dt_local = ds ./ sp;

        % cumulative time (M+1 x 1)
        cumt = [0; cumsum(dt_local)];

        % rescale to actual period of this curve
        Tk = periods(k);
        cumt = cumt * (Tk / cumt(end));

        paramGrid{k} = cumt;
        posGrid{k}   = Cw;
    end

    % --- Mean curve: same idea but using mean speed + mean period ---
    dCm  = diff(posMean,1,1);
    ds_m = sqrt(sum(dCm.^2,2));

    sp_m = meanSpeed;            % Mx1 (already computed)
    sp_m(sp_m < 1e-12) = 1e-12;

    dt_m = ds_m ./ sp_m;
    cumt_mean = [0; cumsum(dt_m)];

    meanPeriod = mean(periods);  % you computed earlier too
    cumt_mean = cumt_mean * (meanPeriod / cumt_mean(end));

    paramMean = cumt_mean;

else
    speed_factor = 0.5;
    % --- Non-dynamics mode: uniform progress in normalized arclength parameter ---
    % Use parameter u in [0,1] for all curves, same speed.
    u = linspace(0,1,M+1).';   % (M+1)x1

    for k = 1:Na
        C = aligned{k};
        posGrid{k}   = [C; C(1,:)];
        paramGrid{k} = u;
    end

    posMean   = [meanCurve; meanCurve(1,:)];
    paramMean = u;
end


%% ================================================================
%   SNAPSHOT EXPORT (6 frames) – put THIS right before the while true loop
%% ================================================================
outDir = 'Figures_vdp/arc_param';
if ~exist(outDir,'dir'); mkdir(outDir); end

nFrames = 9;

% One full cycle in the SAME units as your unified parametrization:
% dynamics=true  -> seconds, use meanPeriod
% dynamics=false -> u in [0,1], use 1
if dynamics
    cycleLen = meanPeriod;
else
    cycleLen = 1;
end

frameTimes = linspace(0, cycleLen, nFrames+1);
frameTimes(end) = [];   % avoid duplicating the start point

for f = 1:numel(frameTimes)
    t_snap = frameTimes(f);

    % --- update all sample trackers ---
    for k = 1:Na
        p = paramGrid{k};      % (M+1)x1
        P = posGrid{k};        % (M+1)x2

        Tk   = p(end);
        tloc = mod(t_snap, Tk);

        px = interp1(p, P(:,1), tloc, 'pchip');
        py = interp1(p, P(:,2), tloc, 'pchip');

        set(hCurve(k),'XData',px,'YData',py);
    end

    % --- update mean tracker ---
    Tm     = paramMean(end);
    tloc_m = mod(t_snap, Tm);

    mx = interp1(paramMean, posMean(:,1), tloc_m, 'pchip');
    my = interp1(paramMean, posMean(:,2), tloc_m, 'pchip');

    set(hMean,'XData',mx,'YData',my);

    drawnow;
    if dynamics
        exportgraphics(gcf, fullfile(outDir, sprintf('fig_frechet_mean_dynamic_frame_%02d.pdf', f)), ...
        'ContentType','vector');
    else
        exportgraphics(gcf, fullfile(outDir, sprintf('fig_frechet_mean_non_dynamic_frame_%02d.pdf', f)), ...
        'ContentType','vector');
    end
end

%% ================================================================
%   ONE UNIFIED ANIMATION LOOP
%% ================================================================
while true
    t_global = t_global + speed_factor * dt_play;

    % update sample curve markers
    for k = 1:Na
        p = paramGrid{k};
        P = posGrid{k};

        T = p(end);                 % end of param domain (period or 1)
        tloc = mod(t_global, T);    % wrap

        px = interp1(p, P(:,1), tloc, 'pchip');
        py = interp1(p, P(:,2), tloc, 'pchip');

        set(hCurve(k),'XData',px,'YData',py);
    end

    % update mean marker
    Tm = paramMean(end);
    tloc_m = mod(t_global, Tm);

    mx = interp1(paramMean, posMean(:,1), tloc_m, 'pchip');
    my = interp1(paramMean, posMean(:,2), tloc_m, 'pchip');

    set(hMean,'XData',mx,'YData',my);

    drawnow limitrate;
    pause(dt_play);
end

end
