function vdp_scaled_cycles_tracking()
% Van der Pol: x'=y - c x^3/3 + x, y' = epsilon(a-x)
% Check if periodic solutions for different c are scaled versions
% Supports time or arc-length parametrization
% Shows moving tracking points along curves

clear; close all; clc;

%% --- Parameters
c_vals = linspace(0.5,10,5);   % values of c to test
epsilon = 0.6;
a = 0;                         
x0 = [0.5;0];                  
Tfinal = 400; dt = 0.005;     
M = 300;                        
residual_tol = 1e-2;            
paramMode = 'arc';  % 'time' or 'arc'

%% --- Storage
cycles = cell(numel(c_vals),1);

%% --- Solve system for each c
for k = 1:numel(c_vals)
    c = c_vals(k);
    vdp = @(t,z)[ z(2) - c*z(1).^3/3 + z(1); epsilon*(a - z(1)) ];
    tspan = 0:dt:Tfinal;
    [t,z] = ode45(vdp, tspan, x0);

    % Take last half of simulation (steady state)
    idx = t >= Tfinal/2;
    Z = z(idx,:);

    % Find one period via peaks
    [~,locs] = findpeaks(Z(:,1),'MinPeakProminence',0.05);
    if numel(locs)<2
        error('Not enough peaks for c=%.3f',c);
    end
    i1 = locs(end-1); i2 = locs(end);
    cycle = Z(i1:i2,:);

    % Resample
    switch paramMode
        case 'time'
            tU = linspace(0,1,M);
            cycle_res = interp1(linspace(0,1,size(cycle,1)), cycle, tU);
        case 'arc'
            cycle_res = resampleArcLengthClosed(cycle,M);
        otherwise
            error('Unknown paramMode');
    end
    cycles{k} = cycle_res;
end

%% --- Compare scaling relative to first cycle
ref = cycles{1};
alpha_vals = zeros(numel(c_vals),1);
residuals = zeros(numel(c_vals),1);
scaled_flag = false(numel(c_vals),1);

for k = 1:numel(c_vals)
    C = cycles{k};
    % Center curves
    ref_centered = ref - mean(ref,1);
    C_centered   = C - mean(C,1);
    % Optimal scaling factor alpha
    alpha = sum(ref_centered(:).*C_centered(:)) / sum(ref_centered(:).^2);
    alpha_vals(k) = alpha;
    % Relative residual
    residuals(k) = norm(C_centered - alpha*ref_centered,'fro') / norm(C_centered,'fro');
    % Determine if approx. scaled
    scaled_flag(k) = residuals(k) < residual_tol;
end

%% --- Display results
fprintf('c\talpha\t\tresidual\tscaled\n');
for k = 1:numel(c_vals)
    fprintf('%.2f\t%.3f\t%.4f\t%s\n', c_vals(k), alpha_vals(k), residuals(k), ...
            tern(scaled_flag(k),'Yes','No'));
end

%% --- Plot cycles with moving tracking points
figure('Color','w'); hold on; axis equal; grid on;
colors = lines(numel(c_vals));
hCurve = gobjects(numel(c_vals),1);
hDot   = gobjects(numel(c_vals),1);

% Plot curves
for k = 1:numel(c_vals)
    C = cycles{k};
    plot(C(:,1), C(:,2), 'Color',colors(k,:),'LineWidth',1.5, 'DisplayName',sprintf('c=%.2f',c_vals(k)));
    % Moving tracking point
    hDot(k) = plot(C(1,1), C(1,2), 'o', 'MarkerFaceColor',colors(k,:), 'Color',colors(k,:),'MarkerSize',8, 'HandleVisibility','off');
end
xlabel('x'); ylabel('y'); title('Van der Pol cycles with tracking points');
legend('Location','best');

%% --- Animation loop
while true
    for i = 1:M
        for k = 1:numel(c_vals)
            C = cycles{k};
            set(hDot(k),'XData',C(i,1),'YData',C(i,2));
        end
        drawnow limitrate;
        pause(0.02);
    end
end

end

%% --- Resample closed curve by arc length
function C = resampleArcLengthClosed(seg, M)
    seg2 = [seg; seg(1,:)];  % close curve
    d = diff(seg2,1,1);
    segLen = sqrt(sum(d.^2,2));
    cumLen = [0; cumsum(segLen)];
    L = cumLen(end);
    s = linspace(0,L,M+1);
    x = interp1(cumLen, seg2(:,1), s);
    y = interp1(cumLen, seg2(:,2), s);
    C = [x(1:end-1).' y(1:end-1).'];
end

%% --- Simple ternary
function out = tern(cond,valTrue,valFalse)
    if cond
        out = valTrue;
    else
        out = valFalse;
    end
end