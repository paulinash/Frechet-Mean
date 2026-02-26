function vdp_fm_time_wwo_period
% Iterative Fréchet (Karcher) mean for Van der Pol periodic curves
% Supports:
%   - time parametrization WITHOUT period handling (Code G)
%   - time parametrization WITH period handling (Code I)

clear; close all; clc;

%% ================================================================
% USER OPTIONS
%% ================================================================

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

% true: all curves are rescaled back to their respective period    
% false: all curves are sampled with M points and thuse scaled onto one
% period
dynamics = false;   
useColor = true;
a = -1.0;
b = 1.0;
Na = 50;
random_samples = true;
nFrames = 0;                  % number of snapshots, 0 if to be skipped


epsilon   = 0.1;

% Samples of a
if random_samples
    rng(123);
    a_vals = a + (b - a) * rand(1, Na);
else    
    a_vals = linspace(a,b,Na);
end

Tfinal = 400;
dt = 0.005;
M = 300;

maxIter = 30;
tol     = 1e-6;

%% ================================================================
% SIMULATE + EXTRACT ONE PERIOD + RESAMPLE
%% ================================================================
curves_res = cell(Na,1);

% NEW (from Code I)
periods = zeros(Na,1);

for k = 1:Na
    a = a_vals(k);
    vdp = @(t,z) [z(2) - z(1)^3/3 + z(1);
                  epsilon*(a - z(1))];

    [t,z] = ode45(vdp, 0:dt:Tfinal, [-1; 0]);

    % remove transient
    idx = t >= Tfinal/2;
    z2 = z(idx,:);
    t2 = t(idx);

    % detect period
    x = z2(:,1);
    [~,locs] = findpeaks(x,'MinPeakProminence',0.05,'MinPeakDistance',10);
    i1 = locs(end-1);
    i2 = locs(end);

    seg  = z2(i1:i2,:);
    tseg = t2(i1:i2);

    periods(k) = tseg(end) - tseg(1);
    fprintf('Period of curve %d: %.2f\n', k, periods(k));


    % ------------------------------------------------------------
    % RESAMPLING
    % ------------------------------------------------------------
    tU = linspace(tseg(1), tseg(end), M);
    curves_res{k} = interp1(tseg, seg, tU);
end

%% ================================================================
% ITERATIVE FRÉCHET / KARCHER MEAN
%% ================================================================
aligned = curves_res;

sq   = linspace(0,1,M+1);
sq   = sq(1:end-1);
tauGrid = linspace(0,1,M);
prevEnergy = NaN;

% Initial mean 
meanCurve = zeros(M,2);
for k = 1:Na
    meanCurve = meanCurve + curves_res{k};
end
meanCurve = meanCurve / Na;

for iter = 1:maxIter
    mean_old = meanCurve;

    energy_k = zeros(Na,1);

    % ------------------------------------------------------------
    % ALIGN CURVES TO CURRENT MEAN
    % ------------------------------------------------------------
    for k = 1:Na
        C = curves_res{k};

        bestScore = Inf;
        bestTau   = 0;

        for tau = tauGrid
            sShift = mod(sq + tau,1);
            Csh = interp1(sq, C, sShift,'pchip');

            score = mean(sum((Csh - meanCurve).^2,2));
            if score < bestScore
                bestScore = score;
                bestTau = tau;
            end
        end

        sShift = mod(sq + bestTau,1);
        aligned{k} = interp1(sq, C, sShift,'pchip');

        % Store current energy for stopping criterion
        energy_k(k) = bestScore;
    end
    % Current energy
    currEnergy = mean(energy_k);
    
    % Update Mean
    meanCurve = mean(cat(3,aligned{:}),3);

    % convergence check
    relEnergyDrop = abs(currEnergy - prevEnergy) / (prevEnergy + eps);
    fprintf('Iter %d: rel energy drop = %.3e\n',iter,relEnergyDrop);

    if relEnergyDrop < tol
        fprintf("Converged. \n")
        break
    end
    prevEnergy = currEnergy;
end

% mean period
meanPeriod = mean(periods);  
fprintf("Mean period = %.6f\n", meanPeriod);


%% ================================================================
% PLOT CURVES
%% ================================================================
figure('Color','w','Position',[100 100 1100 700]); hold on; axis equal; grid on;
hold on; axis equal; grid on;

titleStr = sprintf('Fréchet mean');
if dynamics
    titleStr = [titleStr ', with period'];
end
%title(titleStr)

if useColor
    colors = lines(Na);            % one color per sample    
    colors = 0.3*colors + 0.7;   % mix with white
else
    gray = [0.75 0.75 0.75];       % light gray
    colors = repmat(gray, Na, 1);  % same gray for all samples
end

for k = 1:Na
    plot([aligned{k}(:,1);  aligned{k}(1,1)], [aligned{k}(:,2);  aligned{k}(1,2)],'Color',colors(k,:),'LineWidth',1.0, 'DisplayName',sprintf('a=%.2f',a_vals(k)));
end
plot([meanCurve(:,1); meanCurve(1,1)], [meanCurve(:,2); meanCurve(1,2)],'k','LineWidth',2,'DisplayName','mean curve');

%xlabel('x'); ylabel('y');
exportgraphics(gcf,'Figures_vdp/time_param/fig_frechet_mean.pdf','ContentType','vector');


%% ================================================================
% TRACKING POINTS
%% ================================================================
hCurve = gobjects(Na,1);
for k = 1:Na
    colors = lines(Na);            % one color per sample
    colors = 0.5*colors + 0.5;
    hCurve(k) = plot(aligned{k}(1,1), aligned{k}(1,2),'o', ...
        'Color',colors(k,:), 'MarkerFaceColor',colors(k,:), ...
        'MarkerSize',8);
end

hMean = plot(meanCurve(1,1), meanCurve(1,2),'ko', ...
             'MarkerFaceColor','k','MarkerSize',9);

%% ================================================================
% SNAPSHOT EXPORT 
%% ================================================================
outDir = 'Figures_vdp/time_param';
if ~exist(outDir,'dir'); mkdir(outDir); end

frameTimes = linspace(0, meanPeriod, nFrames+1);
frameTimes(end) = [];                          % drop endpoint (same as 0)

for f = 1:numel(frameTimes)
    t_global = frameTimes(f);

    % --- update sample trackers ---
    for k = 1:Na
        if dynamics
            Tk = periods(k);
            s  = mod(t_global, Tk) / Tk;
        else
            s  = mod(t_global, meanPeriod) / meanPeriod;
        end

        idxf = s*(M-1)+1;
        i1 = floor(idxf);
        i2 = mod(i1,M)+1;
        a  = idxf - i1;

        Ci = (1-a)*aligned{k}(i1,:) + a*aligned{k}(i2,:);
        set(hCurve(k),'XData',Ci(1),'YData',Ci(2));
    end

    % --- update mean tracker ---
    s = mod(t_global, meanPeriod) / meanPeriod;
    idxf = s*(M-1)+1;
    i1 = floor(idxf);
    i2 = mod(i1,M)+1;
    a  = idxf - i1;

    Cm = (1-a)*meanCurve(i1,:) + a*meanCurve(i2,:);
    set(hMean,'XData',Cm(1),'YData',Cm(2));

    drawnow;
    exportgraphics(gcf, fullfile(outDir, sprintf('fig_frechet_mean_frame_%02d.pdf', f)), ...
        'ContentType','vector');
end


%% ================================================================
% ANIMATION LOOP
%% ================================================================
t_global = 0;
dt_play  = 0.02;
speed    = 5;

while true
    t_global = t_global + speed*dt_play;

    for k = 1:Na
        if dynamics
            Tk = periods(k);
            s  = mod(t_global, Tk) / Tk;
        else
            % abstract time scale for all curves
            s = mod(t_global, meanPeriod) / meanPeriod;
        end

        idxf = s*(M-1)+1;
        i1 = floor(idxf);
        i2 = mod(i1,M)+1;
        a  = idxf - i1;

        Ci = (1-a)*aligned{k}(i1,:) + a*aligned{k}(i2,:);
        set(hCurve(k),'XData',Ci(1),'YData',Ci(2));
    end

    % mean tracker
    s = mod(t_global,meanPeriod)/meanPeriod;

    idxf = s*(M-1)+1;
    i1 = floor(idxf);
    i2 = mod(i1,M)+1;
    a  = idxf - i1;
    Cm = (1-a)*meanCurve(i1,:) + a*meanCurve(i2,:);
    set(hMean,'XData',Cm(1),'YData',Cm(2));

    drawnow limitrate;
    pause(dt_play);
end
end



