function ml_frechet_mean_bursting_3D
% ================================================================
% Simulates Morris-Lecar bursting for multiple applied currents.
% Extracts the periods from each simulation using the slow variable.
% Resamples the burst trajectory along normalized arc length.
% Computes the geometrical Fréchet mean of all curves, allowing rotational alignment.
% Reconstructs time parametrization for dynamics
% Plots individual and mean bursts in 2D and 3D, with an animated tracker along the curves.
% ================================================================

clear; close all; clc;

%% ================= USER SETTINGS =================
low = 45;
high = 45.2;
N = 15;
I_vals = (high-low).*rand(N,1) + low;
I_vals = sort(I_vals);
I_vals = linspace(45,45.2,N);
useColor = true; 

Tfinal = 8000;
dt_sim = 0.05;
tspan = 0:dt_sim:Tfinal;

% M needs to be quite high so that the phase shift actually does something
% (>1000)
M = 2000;                        % Number of points along burst curves
maxIter = 40;                    % Max iterations for Fréchet mean
tol = 1e-8;                      % Convergence tolerance

%% ================= MODEL PARAMETERS ================
Vk  = -84; Vl = -60; Vca = 120;
V1  = -1.2; V2 = 18; V3 = 12; V4 = 17.4;
gk  = 8; gl = 2; c = 20; gca = 4; gkca = 0.25;
y0  = 10; phi = 0.23; mu = 0.2; epsilon = 0.005;

basep = struct('Vk',Vk,'Vl',Vl,'Vca',Vca,'V1',V1,'V2',V2,'V3',V3,'V4',V4,...
               'gk',gk,'gl',gl,'c',c,'gca',gca,'gkca',gkca,'y0',y0,...
               'phi',phi,'mu',mu,'eps',epsilon);

%% ================= STORAGE =================
curves = cell(N,1);        % Arc-length sampled curves
vels   = cell(N,1);        % Corresponding velocities
periods = zeros(N,1);      
rawT   = cell(N,1);        % Raw segment times
rawZ   = cell(N,1);        % Raw segment states
spikeCounts = zeros(N,1);  % Number of spikes per burst

%% ================= SIMULATE & EXTRACT BURST PERIOD =================
for k = 1:N
    p = basep; p.I = I_vals(k);

    z0 = [-20; 0; 6];
    [t,Z] = ode15s(@(t,z) ml_rhs(t,z,p), tspan, z0);
    V = Z(:,1); w = Z(:,2); y = Z(:,3);

    % ---- remove transient ----
    idx = t > Tfinal/2;
    t2 = t(idx); V2 = V(idx); w2 = w(idx); y2 = y(idx);

    %% ================= Period detection on low decline =================
    % Use slow variable for peak detection
    y_min = min(y2); y_max = max(y2);      
    yth = y_min + 0.15*(y_max - y_min);    % Low threshold on decline
    dy_s = gradient(y2, t2);               % To detect downward crossings

    % Candidates: downward crossings of threshold
    downIdx = find((y2(1:end-1) > yth) & (y2(2:end) <= yth)) + 1; 
    % Choose indices with actually decreasing gradient
    downIdx = downIdx(downIdx <= numel(dy_s) & dy_s(downIdx) < -1e-4);

    if numel(downIdx) < 2
        error('Not enough valid downward crossings for I=%.2f.', p.I);
    end

    % Extract last burst period strictly
    i1 = downIdx(end-1); i2 = downIdx(end);
    segT = t2(i1:i2);
    segZ = [V2(i1:i2), w2(i1:i2), y2(i1:i2)];

    rawT{k} = segT;
    rawZ{k} = segZ;
    periods(k) = segT(end) - segT(1);

    %% ================= Voltage spike detection =================
    [Vpk, locs, ~, proms] = findpeaks(segZ(:,1),'MinPeakProminence',0.2);
    if isempty(Vpk)
        vlocs = []; spikeCounts(k)=0;
    else
        prom_thr = max(0.25*max(proms),0.2); % consider only prominent spikes
        keep = proms >= prom_thr;
        vlocs = locs(keep);
        spikeCounts(k) = numel(vlocs);
    end
    fprintf('I=%.2f : burst dur=%.3f s, spikes=%d\n', p.I, periods(k), spikeCounts(k));

    %% ================= Arc-length resampling: curve and velocity =================
    % Velocities
    Vdot = gradient(segZ(:,1), segT);
    Wdot = gradient(segZ(:,2), segT);
    Ydot = gradient(segZ(:,3), segT);
    Vraw = [Vdot, Wdot, Ydot];

    dZ = diff(segZ,1,1);
    s = [0; cumsum(sqrt(sum(dZ.^2,2)))];  % physical arc length
    s = s / s(end);                      % normalize to [0,1]
    sq = linspace(0,1,M+1);              % common domain for all curves


    % Interpolate curves and velocities from time to uniform sampled arclength
    curves{k} = [interp1(s, segZ(:,1), sq(1:end-1),'pchip')', ...
                 interp1(s, segZ(:,2), sq(1:end-1),'pchip')', ...
                 interp1(s, segZ(:,3), sq(1:end-1),'pchip')'];
    vels{k}   = interp1(s, Vraw, sq(1:end-1),'pchip');
end

% Validate spike count
if numel(unique(spikeCounts)) ~= 1
    error('ABORT: bursts have different spike counts: %s', mat2str(spikeCounts));
end
fprintf('All bursts validated: spikes=%d\n', spikeCounts(1));
fprintf('Mean Period = %.2f \n', mean(periods));

%% ================= GEOMETRICAL FRECHET MEAN =================
% Using arc-length parametrization and rotational shift
meanC = mean(cat(3,curves{:}),3);  % Initial mean
%meanC = zeros(M,3); % this initial mean sometimes leads to problems in shift
%meanC = [sin(2*pi*(0:M-1)'/M), cos(2*pi*(0:M-1)'/M), sin(4*pi*(0:M-1)'/M)];
aligned = curves;                  % Will store aligned curves
alignedVel = vels;                 % Aligned velocities
tauGrid = linspace(0,1,M);           % rotations in [0,1]

for iter = 1:maxIter % Iterate until convergence
    mean_old = meanC;
    for k = 1:N % Iterate over each curve
        % Search for best alignment of curve k to current mean 
        % using optimal tau
        bestErr = Inf; bestTau = 0;
        for tau = tauGrid % Iterative over tau
            % Rotate curve along arc length 
            sq_shifted = mod(sq + tau, 1)';
            sq_shifted = sq_shifted(1:end-1);  
            % Interpolate shifted curve to arclength grid
            gammaShift = interp1(sq(1:end-1)', curves{k}, sq_shifted, 'pchip', 'extrap');
            % Pointwise distance to current mean
            err = mean(sum((gammaShift - meanC).^2, 2)); 
            if err < bestErr
                % Choose best alignment
                bestErr = err;
                bestTau = tau;
                bestCurve = gammaShift;
            end
        end
        % Fix optimal aligned curve and associated velocity
        aligned{k} = bestCurve;
        sq_shifted_best = mod(sq + bestTau, 1)';
        sq_shifted_best = sq_shifted_best(1:end-1);
        alignedVel{k} = interp1(sq(1:end-1), vels{k}, sq_shifted_best, 'pchip', 'extrap');
        fprintf('Curve %d best rotation tau = %.4f (error = %.6f)\n', k, bestTau, bestErr);
    end
    % Average optimal aligned curves
    meanC = mean(cat(3,aligned{:}),3);

    % Check convergence
    if norm(meanC(:)-mean_old(:)) / (norm(mean_old(:))+eps) < tol
        break
    end
end
fprintf('Fréchet mean converged in %d iterations\n', iter);

%% ================= TIME RECONSTRUCTION FOR ANIMATION =================
cumTime = cell(N,1); posInterp = cell(N,1); allSpeeds = zeros(N,M);

for k=1:N
    C = aligned{k}; Vv = alignedVel{k};
    Cw = [C; C(1,:)];               % Wrap curve for periodicity
    dC = diff(Cw,1,1); ds = sqrt(sum(dC.^2,2)); % Arc-length increments
    sp = sqrt(sum(Vv.^2,2)); sp(sp<1e-12)=1e-12; % Speed magnitudes
    allSpeeds(k,:) = sp.';
    dt_local = ds ./ sp; % Local time increments (time=distance/speed)
    cumt = [0;cumsum(dt_local)]; % Cumulative time
    cumt = cumt * (periods(k)/cumt(end)); % Rescale to curve specific burst duration
    cumTime{k} = cumt; posInterp{k} = Cw;
end
% Same for mean curve
meanSp = mean(allSpeeds,1).'; % Mean speed vector
Cm_w = [meanC; meanC(1,:)];         % Wrap mean curve
dC = diff(Cm_w,1,1); ds = sqrt(sum(dC.^2,2));
spm = meanSp; spm(spm<1e-12)=1e-12;
dt_mean = ds ./ spm;
cumTime_mean = [0;cumsum(dt_mean)];
cumTime_mean = cumTime_mean * mean(periods) / cumTime_mean(end);

if useColor
    colors = lines(N);            % one color per sample
    colors = 0.5*colors + 0.5;   % mix with white

else
    gray = [0.75 0.75 0.75];       % light gray
    colors = repmat(gray, N, 1);  % same gray for all samples
end

%% ================= 2D PLOTS =================
% V(t)
figure('Color','w','Position',[100 100 1400 600]); hold on;
for k=1:N 
    plot(rawT{k} - rawT{k}(1),rawZ{k}(:,1), 'DisplayName', sprintf('I=%.2f', I_vals(k)),'Color',colors(k,:)); 
end
plot(cumTime_mean(1:end-1), Cm_w(1:end-1,1), 'LineWidth',1.5, 'Color','k', 'DisplayName','mean');
xlabel('t'); ylabel('V'); title('V(t) over one burst');
exportgraphics(gcf,'Figures_ml/2D_V.pdf','ContentType','vector');

% y(t)
figure('Color','w','Position',[100 100 1400 600]); hold on;
for k=1:N
    plot(rawT{k}-rawT{k}(1),rawZ{k}(:,3), 'DisplayName',sprintf('I=%.2f', I_vals(k)),'Color',colors(k,:)); 
end
plot(cumTime_mean(1:end-1), Cm_w(1:end-1,3), 'LineWidth',1.5, 'Color','k','DisplayName','mean')
xlabel('t'); ylabel('y'); title('y(t) over one burst'); 
exportgraphics(gcf,'Figures_ml/2D_y.pdf','ContentType','vector');

% w(t)
figure('Color','w','Position',[100 100 1400 600]); hold on;
for k=1:N
    plot(rawT{k}-rawT{k}(1),rawZ{k}(:,2), 'DisplayName',sprintf('I=%.2f', I_vals(k)),'Color',colors(k,:)); 
end
plot(cumTime_mean(1:end-1), Cm_w(1:end-1,2), 'LineWidth',1.5, 'Color','k','DisplayName','mean')
xlabel('t'); ylabel('y'); title('w(t) over one burst'); 
exportgraphics(gcf,'Figures_ml/2D_w.pdf','ContentType','vector');


%% ================= 3D PLOT WITH TRACKERS =================
figure('Color','w','Position',[200 100 900 700]); hold on; grid on;

for k=1:N
    h = plot3([aligned{k}(:,1); aligned{k}(1,1)],[aligned{k}(:,2); aligned{k}(1,2)],...
        [aligned{k}(:,3); aligned{k}(1,3)],'Color',colors(k,:),...
        'LineWidth',1.2,'DisplayName',sprintf('I=%.2f', I_vals(k))); 
    h.Color(4) = 0.3;
end
plot3(meanC(:,1),meanC(:,2),meanC(:,3),'k','LineWidth',3, 'DisplayName','mean curve');


% Initialize moving trackers
hCurve = gobjects(N,1);
for k=1:N
    hCurve(k) = plot3(aligned{k}(1,1),aligned{k}(1,2),aligned{k}(1,3),'o', ...
        'Color',colors(k,:),'MarkerFaceColor',colors(k,:),'MarkerSize',7, 'HandleVisibility','off');
end
hMean = plot3(meanC(1,1),meanC(1,2),meanC(1,3),'ok',...
    'MarkerSize',10, 'MarkerFaceColor','k','DisplayName','mean tracker');

xlabel('V'); ylabel('w'); zlabel('y'); 
title('3D trajectories with mean tracker'); 
view([-30 20]);
exportgraphics(gcf,'Figures_ml/3D.pdf','ContentType','image','Resolution',600);


%% ================================================================
% SNAPSHOT EXPORT (3D tracker positions at selected times)
%% ================================================================
outDir = 'Figures_ml/';
if ~exist(outDir,'dir'); mkdir(outDir); end

nFrames = 9;  % number of snapshots
meanPeriod = mean(periods);

% Choose snapshot times over one mean period (exclude endpoint so 0 != end)
frameTimes = linspace(0, meanPeriod, nFrames+1);
frameTimes(end) = [];

% (Optional) nicer consistent camera for all snapshots
% view([-30 20]);  % you already set view; keep if you want

for f = 1:numel(frameTimes)
    t_global = frameTimes(f);

    % --- update sample trackers (same logic as animation) ---
    for k = 1:N
        tloc = mod(t_global, cumTime{k}(end));  % each curve's reconstructed period
        px = interp1(cumTime{k}, posInterp{k}(:,1), tloc, 'pchip');
        py = interp1(cumTime{k}, posInterp{k}(:,2), tloc, 'pchip');
        pz = interp1(cumTime{k}, posInterp{k}(:,3), tloc, 'pchip');
        set(hCurve(k),'XData',px,'YData',py,'ZData',pz);
    end

    % --- update mean tracker ---
    tloc = mod(t_global, cumTime_mean(end));
    px = interp1(cumTime_mean, Cm_w(:,1), tloc, 'pchip');
    py = interp1(cumTime_mean, Cm_w(:,2), tloc, 'pchip');
    pz = interp1(cumTime_mean, Cm_w(:,3), tloc, 'pchip');
    set(hMean,'XData',px,'YData',py,'ZData',pz);

    drawnow;

    % Export snapshot
    exportgraphics(gcf, fullfile(outDir, sprintf('3D_frame_%02d.pdf', f)), ...
        'Resolution',600,'ContentType','image');
end



%% ================= ANIMATION LOOP =================
fprintf('Animation running... Press Ctrl+C to stop.\n');
t_global=0; dt_play=0.02; speed_factor=20;

while true
    t_global = t_global + speed_factor*dt_play;

    % Sample curves
    for k=1:N
        tloc = mod(t_global, cumTime{k}(end));
        px = interp1(cumTime{k}, posInterp{k}(:,1), tloc,'pchip');
        py = interp1(cumTime{k}, posInterp{k}(:,2), tloc,'pchip');
        pz = interp1(cumTime{k}, posInterp{k}(:,3), tloc,'pchip');
        set(hCurve(k),'XData',px,'YData',py,'ZData',pz);
    end

    % Sample mean curve
    tloc = mod(t_global, cumTime_mean(end));
    px = interp1(cumTime_mean,Cm_w(:,1),tloc,'pchip');
    py = interp1(cumTime_mean,Cm_w(:,2),tloc,'pchip');
    pz = interp1(cumTime_mean,Cm_w(:,3),tloc,'pchip');
    set(hMean,'XData',px,'YData',py,'ZData',pz);

    drawnow limitrate; pause(dt_play);
end
end

%% ================== RHS ==================
function dz = ml_rhs(~,z,p)
V=z(1); w=z(2); y=z(3);
m_inf=0.5*(1+tanh((V-p.V1)/p.V2));
w_inf=0.5*(1+tanh((V-p.V3)/p.V4));
tau_w=cosh((V-p.V3)/(2*p.V4));
Ica=p.gca*m_inf*(V-p.Vca);
Kfac=p.gk*w + p.gkca*(y/(y+p.y0));
IK=Kfac*(V-p.Vk); Ileak=p.gl*(V-p.Vl);
dV=(1/p.c)*(p.I-Ica-IK-Ileak);
dw=p.phi*tau_w*(w_inf-w);
dy=p.eps*(-p.mu*Ica-y);
dz=[dV; dw; dy];
end
