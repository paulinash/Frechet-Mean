function vdp_frechet_mean_arclength_dynamic_general
% Iterative Fréchet (Karcher) mean for periodic solutions of a **general**
% Van-der-Pol-type system:
%
%   x' = b*y - c*x^3/3 + d*x + e
%   y' = epsilon*(a - x)
%
% You may choose ONE parameter to sample (e.g., vary b from 0 to 2),
% while the remaining parameters stay fixed in their baseline values:
%
%   a = 0,  b = 1,  c = 1,  d = 1,  e = 0
%
% This script then computes periodic solutions, aligns them by arclength,
% computes a Fréchet mean curve, and visualizes with a velocity-driven animation.

clear; close all; clc;

%% ===================== USER OPTIONS ===============================
epsilon = 0.1;    
% Frechet mean iteration settings
maxIter = 40;
tol     = 1e-8;

% ---------- BASELINE VALUES ----------
a_base = 0;
b_base = 1;
c_base = 1;
d_base = 1;
e_base = 0;

% ---------- PICK ONE PARAMETER TO SAMPLE ----------
%   Options: 'a', 'b', 'c', 'd', 'e'
paramToVary = 'd';

switch paramToVary
    case 'a'
        param_vals = linspace(-0.8,0.5,7);
    case 'b' % b in (0.2,10) sensible, changes characteristics at b=1
        param_vals = linspace(0.2,1,10);
    case 'c'
        param_vals = linspace(0.1,3,5);
    case 'd' % sensible in (0.1,5)
        param_vals = linspace(0.1,2,5);
    case 'e'
        param_vals = linspace(-1,1,5);
    otherwise
        error('Invalid parameter name.');
end

Na = numel(param_vals);

%% Simulation time and discretization
Tfinal = 400;
dt = 0.005;
M = 300;

%% STORAGE
curves_res = cell(Na,1);
periods = zeros(Na,1);
vel_time = cell(Na,1);
raw_time_samples = cell(Na,1);

%% ================================================================
%   SIMULATE FOR EACH PARAMETER SAMPLE
%% ================================================================
for k = 1:Na

    % ----- ASSIGN PARAMETERS -----
    a = a_base; b = b_base; c = c_base; d = d_base; e = e_base;

    switch paramToVary
        case 'a', a = param_vals(k);
        case 'b', b = param_vals(k);
        case 'c', c = param_vals(k);
        case 'd', d = param_vals(k);
        case 'e', e = param_vals(k);
    end

    % ----- GENERALIZED VAN-DER-POL ODE -----
    vdp = @(t,z) [ ...
        b*z(2) - c*z(1)^3/3 + d*z(1) + e ; ...
        epsilon*(a - z(1)) ...
    ];
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
    cumLen = [0; cumsum(segLen)]; L = cumLen(end);

    % interpolate curves and velocities to match M points and arc length
    s = linspace(0, L, M+1);
    x_s = interp1(cumLen, seg(:,1), s);
    y_s = interp1(cumLen, seg(:,2), s);
    C = [x_s(1:end-1).' y_s(1:end-1).'];

    vx_s = interp1(cumLen, vx, s);
    vy_s = interp1(cumLen, vy, s);
    V_s = [vx_s(1:end-1).' vy_s(1:end-1).'];

    curves_res{k} = C;
    vel_time{k} = V_s;
end

%% ================================================================
%   ITERATIVE FRÉCHET / KARCHER MEAN
%% ================================================================
aligned = curves_res;
alignedVel = vel_time;

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
        C = curves_res{k}; bestScore = Inf; bestShift = 0;

        % Find best alignment to mean curve for each curve
        for sft = 0:M-1
            Csh = circshift(C,-sft);
            score = mean(sum((Csh-meanCurve).^2,2));
            if score < bestScore
                bestScore = score; 
                bestShift = sft; 
            end
        end
        aligned{k}    = circshift(C, -bestShift);
        alignedVel{k} = circshift(vel_time{k}, -bestShift);
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
fprintf("Mean period = %.6f\n", meanPeriod);

% Compute raw mean velocity
meanVel_raw = zeros(M,2);
for k = 1:Na
    meanVel_raw = meanVel_raw + alignedVel{k}; 
end
meanVel_raw = meanVel_raw / Na;

% Compute tangent vector of mean curve 
tx = gradient(meanCurve(:,1)); ty = gradient(meanCurve(:,2));
T = [tx ty]; Tnorm = sqrt(sum(T.^2,2))+eps;
unitT = T./Tnorm;

% Projecit raw mean velocty onto tangent direction -> mean tangential speed
meanSpeed_along = sum(meanVel_raw .* unitT, 2);
% Rescale new mean velocity to reobtain original speed
meanVel_geom = unitT .* meanSpeed_along;

%% ================================================================
%   PLOT GEOMETRY + MOVING VELOCITY MARKERS
%% ================================================================
% Initialize plot
figure('Color','w','Position',[100 100 1100 700]); hold on; axis equal; grid on;
title('Fréchet Mean of periodic Van Der Pol solutions');
colors = lines(Na);

% Plot sample and mean curves
for k = 1:Na
    plot(aligned{k}(:,1), aligned{k}(:,2),'Color',colors(k,:),'LineWidth',1.0, 'DisplayName',sprintf('parameter=%.2f',param_vals(k)));
end
plot(meanCurve(:,1), meanCurve(:,2),'k','LineWidth',2,'DisplayName','mean curve');

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
legend('Location','southeast');

%% ================================================================
%   ANIMATION LOOP (VELOCITY-DRIVEN TRAJECTORIES)
%% ================================================================

% --- Precompute time (cumulative) for each curve so motion follows speed ---
cumTime  = cell(Na,1);
posInterp = cell(Na,1);

% Retransform arc-length parameterization into true time parameterization
% for sample curves
for k = 1:Na
    C = aligned{k};      % M×2
    V = alignedVel{k};   % M×2 

    % wrap curve (M+1×2) so that it is periodic (append first point)
    Cw = [C; C(1,:)];

    % arc-length increments (M×1)
    dC = diff(Cw,1,1);
    ds = sqrt(sum(dC.^2,2));

    % speed magnitude (M×1)
    sp = sqrt(sum(V.^2,2));  
    sp(sp < 1e-12) = 1e-12;

    % local dt local time step (distance/speed) (M×1)
    dt_local = ds ./ sp;

    % cumulative time (M+1×1)
    cumt = [0; cumsum(dt_local)];

    % rescale to actual period
    Tk = periods(k);
    cumt = cumt * (Tk / cumt(end));

    cumTime{k}  = cumt;   % (M+1×1)
    posInterp{k} = Cw;    % (M+1×2)
end

% --- Same for mean curve ---
Cw_mean = [meanCurve; meanCurve(1,:)];
V_mean  = meanVel_geom;   % M×2

dCm = diff(Cw_mean,1,1);
ds_m = sqrt(sum(dCm.^2,2));

sp_m = sqrt(sum(V_mean.^2,2));
sp_m(sp_m < 1e-12) = 1e-12;

dt_m = ds_m ./ sp_m;
cumt_mean = [0; cumsum(dt_m)];
cumt_mean = cumt_mean * (meanPeriod / cumt_mean(end));

fprintf("Animation running... Press Ctrl+C to stop.\n");

%% Actual animation
% initialize time variables for animation
t_global = 0;
dt_play = 0.01;
speed_factor = 5; % choose the animation speed

% Endless animation loop
while true
    % Advancing time variable
    t_global = t_global + speed_factor*dt_play;

    % Update all sample curve markers
    for k = 1:Na
        cumt = cumTime{k};    % M+1×1
        Cw   = posInterp{k};  % M+1×2
        Tk   = cumt(end);

        % wrap local time into [0,Tk) to loop motion
        tloc = mod(t_global, Tk);

        % interpolate position
        px = interp1(cumt, Cw(:,1), tloc, 'pchip');
        py = interp1(cumt, Cw(:,2), tloc, 'pchip');

        % move marker to neew position
        set(hCurve(k),'XData',px,'YData',py);
    end

    % Same procedure for mean curve
    tloc_m = mod(t_global, cumt_mean(end));
    mx = interp1(cumt_mean, Cw_mean(:,1), tloc_m, 'pchip');
    my = interp1(cumt_mean, Cw_mean(:,2), tloc_m, 'pchip');
    set(hMean,'XData',mx,'YData',my);

    drawnow limitrate;
    pause(dt_play);
end
end
