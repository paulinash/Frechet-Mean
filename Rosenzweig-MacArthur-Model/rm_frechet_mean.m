function rm_frechet_mean
% rm_frechet_mean_arclength
% Compute iterative Fréchet mean (arclength parametrization) for several
% Rosenzweig-MacArthur limit cycles and animate moving markers on sample
% curves and the mean curve.
%
% Model (nondimensional):
%   dx/dt = x*(1 - x/gamma) - x*y/(1 + x)
%   dy/dt = beta*( x/(1+x) - alpha )*y
%
% Hopf threshold: gamma_H = (1+alpha)/(1-alpha). Choose gamma > gamma_H.

clear; close all; clc;

% Plotting conventions
set(groot, ...
    'DefaultTextInterpreter','latex', ...
    'DefaultAxesTickLabelInterpreter','latex', ...
    'DefaultLegendInterpreter','latex', ...
    'DefaultAxesFontName','Times', ...
    'DefaultTextFontName','Times', ...
    'DefaultAxesFontSize', 12, ...
    'DefaultTextFontSize', 12, ...
    'DefaultLegendFontSize', 11, ...
    'DefaultLineLineWidth', 1.2, ...
    'DefaultFigureColor','w');


%% ========================== USER OPTIONS =============================
% Choose alpha and beta values to explore. Gamma will be chosen relative
% to the Hopf threshold gamma_H to ensure oscillations.

alphaN = 10;
alpha_low = 0.2;
alpha_high = 0.6;

betaN = 1;
beta_low = 0.1;
beta_high = 2;

gamma_multN = 1;
gamma_mult_low = 1.1;
gamma_mult_high = 2;

% Random seed
rng(123);

alpha_vals = alpha_low + (alpha_high - alpha_low)*rand(1,alphaN);
beta_vals = beta_low + (beta_high - beta_low)*rand(1,betaN);
gamma_mults = gamma_mult_low + (gamma_mult_high - gamma_mult_low)*rand(1,gamma_multN);


% Compose parameter list (all combinations)
paramList = [];
for ia = 1:numel(alpha_vals)
    for ib = 1:numel(beta_vals)
        for ig = 1:numel(gamma_mults)
            alpha = alpha_vals(ia);
            beta  = beta_vals(ib);
            gamma = ((1+alpha)/(1-alpha)) * gamma_mults(ig); % > gamma_H
            paramList(end+1).alpha = alpha; 
            paramList(end).beta  = beta;
            paramList(end).gamma = gamma;
        end
    end
end

Na = numel(paramList);
fprintf("Will simulate %d parameter cases.\n", Na);
for k=1:Na
    alpha = paramList(k).alpha;
    gamma = paramList(k).gamma;
    gammaH = (1+alpha)/(1-alpha);
    fprintf("Case %2d: alpha=%.3f, beta=%.3f, gamma=%.3f (gamma_H=%.3f, mult=%.2f)\n", ...
        k, paramList(k).alpha, paramList(k).beta, paramList(k).gamma, gammaH, gamma/gammaH);
end

% Simulation time / solver settings
Tfinal = 2000;      % make long enough to let transients die
dt     = 0.02;
tspan  = 0:dt:Tfinal;
M      = 600;      % number of geometric samples per curve (perimeter points)

% Frechet mean settings
maxIter = 20;
tol     = 1e-8;

% Snapshot settings
nFrames = 6;


% Peak detection settings (to extract last single period)
minProm = 1e-3;
minPeakDist_seconds = 5;   % minimal separation between peaks in time units
minPeakDist_pts = round(minPeakDist_seconds / dt);

%% STORAGE
curves_res = cell(Na,1);    % M x 2 arclength-sampled curves
vel_time   = cell(Na,1);    % M x 2 velocities (time derivatives) sampled at those points
periods    = zeros(Na,1);   % period per case
raw_time_samples = cell(Na,1);

%% ====================== SIMULATE + EXTRACT PERIOD ======================
for k = 1:Na
    alpha = paramList(k).alpha;
    beta  = paramList(k).beta;
    gamma = paramList(k).gamma;

    % Define ODE
    rm = @(t,z) [
        z(1).*(1 - z(1)/gamma) - z(1).*z(2)./(1+z(1));
        beta.*( z(1)./(1+z(1)) - alpha ).*z(2)
        ];

    % initial condition (small predator & prey positive)
    z0 = [0.5;0.5];
    % Integrate
    [t,z] = ode45(rm, tspan, z0);

    % discard first half (transients)
    idx = t >= (Tfinal/2);
    t2 = t(idx); z2 = z(idx,:);
    x = z2(:,1);

    % Find peaks in prey to determine period (take last two peaks)
    [pks,locs] = findpeaks(x,'MinPeakProminence',minProm,'MinPeakDistance',minPeakDist_pts);
    if numel(locs) < 2
        error("Not enough peaks for case %d (alpha=%.3f,beta=%.3f,gamma=%.3f). Adjust time or parameters.", k, alpha, beta, gamma);
    end
    % use last two peaks from the retained segment
    i1 = locs(end-1); i2 = locs(end);

    seg  = z2(i1:i2,:);       % one full period (positions)
    tseg = t2(i1:i2);        % corresponding times (not uniformly spaced)
    raw_time_samples{k}.pos = seg;
    raw_time_samples{k}.t   = tseg;
    Tk = tseg(end) - tseg(1);
    periods(k) = Tk;
    fprintf('Curve %d: period %.2f \n', k, periods(k));

    % Compute velocities (time derivatives) along the extracted segment
    vx = gradient(seg(:,1), tseg);
    vy = gradient(seg(:,2), tseg);

    % compute cumulative arc-length along the segment
    d = diff(seg,1,1); segLen = sqrt(sum(d.^2,2));
    cumLen = [0; cumsum(segLen)];
    s_norm = cumLen / cumLen(end);    % normalize to [0,1]

    sq = linspace(0,1,M+1);           % common domain

    x_s = interp1(s_norm, seg(:,1), sq, 'pchip');
    y_s = interp1(s_norm, seg(:,2), sq, 'pchip');

    C = [x_s(1:end-1).' y_s(1:end-1).'];   % M x 2 (drop duplicate endpoint)

    vx_s = interp1(s_norm, vx, sq, 'pchip');
    vy_s = interp1(s_norm, vy, sq, 'pchip');
    V_s = [vx_s(1:end-1).' vy_s(1:end-1).'];

    curves_res{k} = C;
    vel_time{k}   = V_s;
end

%% ====================== ITERATIVE FRÉCHET / KARCHER MEAN ======================
aligned = curves_res;
alignedVel = vel_time;
tauGrid = linspace(0,1,M);

% Initial guess: simple pointwise average (on the arclength grid)
meanCurve = zeros(M,2);
for k = 1:Na
    meanCurve = meanCurve + curves_res{k};
end
meanCurve = meanCurve / Na;

% === ADDED: Initial alignment to mu_0 to get true F^(N)(mu_0) ===
energy_k = zeros(Na,1);
for k = 1:Na
    C = curves_res{k}; bestScore = Inf;
    for tau = tauGrid
        sq_shift = mod(sq(1:end-1) + tau, 1);
        Csh = interp1(sq(1:end-1), C, sq_shift, 'pchip');
        score = mean(sum((Csh - meanCurve).^2, 2));
        if score < bestScore
            bestScore = score;
            bestCurve = Csh;
            bestTau = tau;
        end
    end
    aligned{k} = bestCurve;
    sq_shift = mod(sq(1:end-1) + bestTau, 1);
    alignedVel{k} = interp1(sq(1:end-1), vel_time{k}, sq_shift, 'pchip');
    energy_k(k) = bestScore;
end
prevEnergy = mean(energy_k);

% Karcher mean iterations: align by circular shift (time/phase shift on the arclength grid)
for iter = 1:maxIter

% === MOVED: Recompute mean from aligned curves first ===
    meanCurve = zeros(M,2);
    for k = 1:Na
        meanCurve = meanCurve + aligned{k};
    end
    meanCurve = meanCurve / Na;

    energy_k = zeros(Na,1);
% Re-align each curve to new mean by best circular shift
    for k = 1:Na
        C = curves_res{k}; bestScore = Inf;
        for tau = tauGrid
            sq_shift = mod(sq(1:end-1) + tau, 1);
            Csh = interp1(sq(1:end-1), C, sq_shift, 'pchip');
            score = mean(sum((Csh - meanCurve).^2, 2));
            if score < bestScore
                bestScore = score;
                bestCurve = Csh;
                bestTau = tau;
            end
        end
        aligned{k} = bestCurve;
        sq_shift = mod(sq(1:end-1) + bestTau, 1);
        alignedVel{k} = interp1(sq(1:end-1), vel_time{k}, sq_shift, 'pchip');
        energy_k(k) = bestScore;
    end
    currEnergy = mean(energy_k);

% convergence check
    relEnergyDrop = abs(currEnergy - prevEnergy) / (prevEnergy+eps);
    fprintf('Iter %d: rel change = %.3e\n', iter, relEnergyDrop);
    if relEnergyDrop < tol
        fprintf('Converged at iter %d.\n', iter);
        break;
    end
    prevEnergy = currEnergy;
end


%% ====================== COMPUTE MEAN PERIOD & MEAN VELOCITY ======================
meanPeriod = mean(periods);
fprintf('Mean period (average over cases) = %.6f\n', meanPeriod);

allSpeeds = zeros(Na,M);
for k = 1:Na
    V = alignedVel{k};
    allSpeeds(k,:) = sqrt(sum(V.^2,2));
end

meanSpeed = mean(allSpeeds,1).';   % M x 1
%meanSpeed = exp(mean(log(allSpeeds),1)).'; geometric mean

%% ====================== PREPARE FOR ANIMATION (TIME-REPARAMETRIZE) ==========
% For each aligned curve, rebuild a periodic (M+1) array and compute
% cumulative time along curve using local ds / speed so animation moves
% along each cycle with physically consistent speed.

cumTime  = cell(Na,1);
posInterp = cell(Na,1);

for k = 1:Na
    C = aligned{k};      % M x 2
    V = alignedVel{k};   % M x 2

    % wrap curve (M+1 x 2)
    Cw = [C; C(1,:)];

    % arc-length increments (M x 1)
    dC = diff(Cw,1,1);
    ds = sqrt(sum(dC.^2,2));

    % speed magnitude from velocity samples (M x 1)
    sp = sqrt(sum(V.^2,2));
    sp(sp < 1e-12) = 1e-12;

    % local dt approximated by ds / speed
    dt_local = ds ./ sp;
    cumt = [0; cumsum(dt_local)];
    % rescale cumt to the measured period for this case
    Tk = periods(k);
    if cumt(end) <= 0
        error('Nonpositive cumt end for case %d', k);
    end
    cumt = cumt * (Tk / cumt(end));
    cumTime{k} = cumt;         % M+1 x 1
    posInterp{k} = Cw;         % M+1 x 2
end

% mean curve wrap
Cw_mean = [meanCurve; meanCurve(1,:)];
dCm = diff(Cw_mean,1,1);
ds_m = sqrt(sum(dCm.^2,2));

sp_m = meanSpeed;
sp_m(sp_m < 1e-12) = 1e-12;
dt_m = ds_m ./ sp_m;
cumt_mean = [0; cumsum(dt_m)];
cumt_mean = cumt_mean * (meanPeriod / cumt_mean(end));
mean_period_intrinsic = sum(ds_m ./ sp_m);
fprintf('Intrinsic mean period: %.4f \n', mean_period_intrinsic);

%% ====================== PLOT + ANIMATE ==================================
figure('Color','w','Position',[50 50 1200 700]); hold on; axis equal; grid on;
colors = lines(Na);
colors = 0.3*colors + 0.7;   % mix with white


% plot sample curves and mean
for k = 1:Na
    plot([aligned{k}(:,1); aligned{k}(1,1)], [aligned{k}(:,2); aligned{k}(1,2)], 'Color', colors(k,:), 'LineWidth', 1.2, ...
        'DisplayName', sprintf('case %d: a=%.2f b=%.2f g=%.2f', k, paramList(k).alpha, paramList(k).beta, paramList(k).gamma)); 
end
plot([meanCurve(:,1); meanCurve(1,1)], [meanCurve(:,2); meanCurve(1,2)], 'k', 'LineWidth', 2.2, 'DisplayName','mean curve');

% File naming
filename = 'Figures_rm/fig_frechet_mean';
if alphaN > 1
    filename = sprintf('%s_alpha', filename);
end
if betaN > 1
    filename = sprintf('%s_beta', filename);
end
if gamma_multN > 1
    filename = sprintf('%s_gamma', filename);
end
filename = sprintf('%s.pdf', filename);

exportgraphics(gcf,filename,'ContentType','vector');

% create markers for animation
hCurve = gobjects(Na,1);
for k = 1:Na
    C = aligned{k};
    colors = lines(Na);
    colors = 0.5*colors + 0.5;   % mix with white
    hCurve(k) = plot(C(1,1), C(1,2), 'o', 'Color', colors(k,:), ...
        'MarkerFaceColor', colors(k,:), 'MarkerSize', 7, 'HandleVisibility','off');
end
hMean = plot(meanCurve(1,1), meanCurve(1,2), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 9, 'DisplayName','tracker');
fprintf('Animation running... press Ctrl-C in the MATLAB window to stop.\n');

%% ================================================================
%   SNAPSHOT EXPORT 
%% ================================================================

frameTimes = linspace(0, meanPeriod, nFrames+1);
frameTimes(end) = [];   % avoid duplicating the start point

for f = 1:numel(frameTimes)
    t_snap = frameTimes(f);

    % --- update all sample trackers ---
    for k = 1:Na
        p = cumTime{k};      % (M+1)x1
        P = posInterp{k};        % (M+1)x2

        Tk   = p(end);
        tloc = mod(t_snap, Tk);

        px = interp1(p, P(:,1), tloc, 'pchip');
        py = interp1(p, P(:,2), tloc, 'pchip');

        set(hCurve(k),'XData',px,'YData',py);
    end

    % --- update mean tracker ---
    Tm     = cumt_mean(end);
    tloc_m = mod(t_snap, Tm);

    mx = interp1(cumt_mean, Cw_mean(:,1), tloc_m, 'pchip');
    my = interp1(cumt_mean, Cw_mean(:,2), tloc_m, 'pchip');

    set(hMean,'XData',mx,'YData',my);

    drawnow;
    % File naming
    filename = 'Figures_rm/fig_frechet_mean';
    if alphaN > 1
        filename = sprintf('%s_alpha', filename);
    end
    if betaN > 1
        filename = sprintf('%s_beta', filename);
    end
    if gamma_multN > 1
        filename = sprintf('%s_gamma', filename);
    end
    filename = sprintf('%s_snapshot_%02d.pdf', filename, f);

    exportgraphics(gcf, filename,'ContentType','vector');
end


%% Endless animation loop
t_global = 0;
dt_play = 0.02;
speed_factor = 1.5;

while true
    t_global = t_global + speed_factor*dt_play;

    for k = 1:Na
        cumt = cumTime{k};
        Cw   = posInterp{k};
        Tk   = cumt(end);
        tloc = mod(t_global, Tk);

        px = interp1(cumt, Cw(:,1), tloc, 'pchip');
        py = interp1(cumt, Cw(:,2), tloc, 'pchip');

        set(hCurve(k), 'XData', px, 'YData', py);
    end

    tloc_m = mod(t_global, cumt_mean(end));
    mx = interp1(cumt_mean, Cw_mean(:,1), tloc_m, 'pchip');
    my = interp1(cumt_mean, Cw_mean(:,2), tloc_m, 'pchip');
    set(hMean, 'XData', mx, 'YData', my);

    drawnow limitrate;
    pause(dt_play);
end

end
