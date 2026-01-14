function ml_frechet_mean_bursting_3D
% ================================================================
% Fréchet mean for Morris-Lecar bursting
% - Period detection on the low decline of y
% - Strict segment boundaries (no extension)
% - Adaptive spike detection
% - Full-time and segment diagnostics
% ================================================================

clear; close all; clc;

%% ================= USER SETTINGS =================

I_vals = [45,45.1,45.2];   % example
N = numel(I_vals);

Tfinal = 8000;
dt_sim = 0.05;
tspan = 0:dt_sim:Tfinal;

M = 1000;
maxIter = 40;
tol = 1e-8;

%% ================= MODEL PARAMETERS ================
Vk  = -84; Vl = -60; Vca = 120;
V1  = -1.2; V2 = 18;
V3  = 12; V4 = 17.4;
gk  = 8; gl = 2; c = 20;
gca = 4; gkca = 0.25; y0 = 10; phi = 0.23;
mu = 0.2; epsilon = 0.005;

basep = struct('Vk',Vk,'Vl',Vl,'Vca',Vca,'V1',V1,'V2',V2,'V3',V3,'V4',V4,...
               'gk',gk,'gl',gl,'c',c,'gca',gca,'gkca',gkca,'y0',y0,...
               'phi',phi,'mu',mu,'eps',epsilon);

%% ================= STORAGE =================
curves = cell(N,1);
vels = cell(N,1);
periods = zeros(N,1);
rawT = cell(N,1); rawZ = cell(N,1);
spikeCounts = zeros(N,1);

%% ================= SIMULATE & EXTRACT BURST PERIOD =================
for k = 1:N
    p = basep;
    p.I = I_vals(k);

    z0 = [-20; 0; 6];
    [t,Z] = ode15s(@(t,z) ml_rhs(t,z,p), tspan, z0);

    V = Z(:,1); w = Z(:,2); y = Z(:,3);

    % ---- full-time plot ----
    figure('Name',sprintf('Full time I=%.2f',p.I),'Color','w');
    subplot(2,1,1); plot(t,V); xlabel('t'); ylabel('V'); title(sprintf('Full V(t), I=%.2f',p.I)); grid on;
    subplot(2,1,2); plot(t,y); xlabel('t'); ylabel('y'); title(sprintf('Full y(t), I=%.2f',p.I)); grid on;

    % remove transient
    idx = t > Tfinal/2;
    t2 = t(idx); V2 = V(idx); w2 = w(idx); y2 = y(idx);

   
    %% ================= Period detection on low decline =================
    y_min = min(y2); y_max = max(y2);
    alpha = 0.15; % low threshold on decline
    yth = y_min + alpha*(y_max - y_min);

    dy_s = gradient(y2, t2);
    % Candidates: y crosses threshold downward (decline)
    downIdx = find( (y2(1:end-1) > yth) & (y2(2:end) <= yth) ) + 1;

    % Ensure derivative negative at crossing (decline)
    downIdx = downIdx( downIdx <= numel(dy_s) & dy_s(downIdx) < -1e-4 );

    if numel(downIdx) < 2
        error('Not enough valid downward crossings for I=%.2f.', p.I);
    end

    % Use last two crossings
    i1 = downIdx(end-1);
    i2 = downIdx(end);

    % Extract segment strictly between period markers (no extension)
    segT = t2(i1:i2);
    segZ = [V2(i1:i2), w2(i1:i2), y2(i1:i2)];

    rawT{k} = segT;
    rawZ{k} = segZ;

    periods(k) = segT(end) - segT(1);

    %% ================= Spike detection =================
    [Vpk, locs, wpk, proms] = findpeaks(segZ(:,1), 'MinPeakProminence', 0.2);
    if isempty(Vpk)
        vlocs = [];
        spikeCounts(k) = 0;
    else
        prom_thr = max(0.25 * max(proms), 0.2);
        keep = proms >= prom_thr;
        vlocs = locs(keep);
        spikeCounts(k) = numel(vlocs);
    end

    fprintf('I=%.2f : burst dur=%.3f s, spikes=%d (raw %d)\n', ...
        p.I, periods(k), spikeCounts(k), numel(locs));

    %% ================= Diagnostic plotting =================
    figure('Name',sprintf('Debug extraction I=%.2f',p.I),'Color','w');
    subplot(3,1,1); hold on;
    plot(t2, y2, 'LineWidth', 1);
    yline(yth,'k--','LineWidth',1.2);
    plot(t2(downIdx), y2(downIdx),'ro','MarkerFaceColor','r');
    plot([segT(1) segT(end)], [min(y2) min(y2)], 'g', 'LineWidth',3);
    xlabel('t'); ylabel('y2'); title('y2 and downward crossings'); grid on;

    subplot(3,1,2); hold on;
    plot(segT, segZ(:,1),'b');
    if ~isempty(locs)
        plot(segT(locs), segZ(locs,1),'x','MarkerSize',8,'LineWidth',1.2,'Color',[0.8 0.2 0.2]);
    end
    if ~isempty(vlocs)
        plot(segT(vlocs), segZ(vlocs,1),'o','MarkerFaceColor','g','MarkerEdgeColor','k');
    end
    xlabel('t'); ylabel('V'); title('Segment V(t) with spikes'); grid on;

    subplot(3,1,3); hold on;
    plot(segT, segZ(:,3));
    plot(segT(1), segZ(1,3),'go','MarkerFaceColor','g');
    plot(segT(end), segZ(end,3),'ro','MarkerFaceColor','r');
    xlabel('t'); ylabel('y'); title('Segment y(t)'); grid on;

    % ====== Arc-length resampling ======
    vdot = gradient(segZ(:,1), segT);
    wdot = gradient(segZ(:,2), segT);
    ydot = gradient(segZ(:,3), segT);
    Vraw = [vdot wdot ydot];

    dZ = diff(segZ,1,1); ds = sqrt(sum(dZ.^2,2));
    s = [0; cumsum(ds)]; L = s(end);
    sq = linspace(0,L,M+1);

    Xq = interp1(s, segZ(:,1), sq, 'pchip');
    Wq = interp1(s, segZ(:,2), sq, 'pchip');
    Yq = interp1(s, segZ(:,3), sq, 'pchip');

    curves{k} = [Xq(1:end-1).' Wq(1:end-1).' Yq(1:end-1).'];
    Vq = interp1(s, Vraw, sq, 'pchip');
    vels{k} = Vq(1:end-1,:);
end

%% ================= Validate spike counts =================
if numel(unique(spikeCounts)) ~= 1
    error('ABORT: bursts have different spike counts: %s', mat2str(spikeCounts));
end
fprintf('All bursts validated: spikes=%d\n', spikeCounts(1));



%% ========== RHS =================
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

end