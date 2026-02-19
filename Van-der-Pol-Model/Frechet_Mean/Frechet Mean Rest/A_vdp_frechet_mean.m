function vdp_frechet_mean_no_phase()
% vdp_frechet_mean_no_phase
% Compute Fréchet mean of Van der Pol limit cycles without phase alignment

clear; close all; clc;

%% Parameters
epsilon = 0.1;                 
a_vals = linspace(-0.9, 0.9, 2); 
Tfinal = 100;                  
dt = 0.01;                      
M = 500;                        
useArcLength = false;  % true = uniform arc-length, false = uniform time

%% Preallocate curves
curves = cell(length(a_vals),1);

%% Solve Van der Pol for each 'a'
for k = 1:length(a_vals)
    a = a_vals(k);
    vdp = @(t,z) [ z(2) - z(1)^3/3 + z(1); epsilon*(a - z(1)) ];
    z0 = [-1;0];
    [t,z] = ode45(vdp,0:dt:Tfinal,z0);
    
    % Take last half of trajectory
    z_ss = z(t >= Tfinal/2,:);
    
    % Resample
    if useArcLength
        curves{k} = resampleArcLength(z_ss, M);
    else
        curves{k} = interp1(linspace(0,1,size(z_ss,1)), z_ss, linspace(0,1,M));
    end
end

fprintf('curve 1: %g\n', curves{1})
fprintf('curve 2: %g\n', curves{2})


%% Compute Fréchet mean (simple pointwise average)
allCurves = zeros(M,2,length(curves));
for k = 1:length(curves)
    allCurves(:,:,k) = curves{k};
end

meanCurve = mean(allCurves,3);

%% Plot curves and mean
figure; hold on; grid on;
colors = lines(length(a_vals));
for k = 1:length(a_vals)
    plot(allCurves(:,1,k), allCurves(:,2,k), 'Color', colors(k,:), 'LineWidth',1.5);
end
plot(meanCurve(:,1), meanCurve(:,2), 'k', 'LineWidth',3);

xlabel('x'); ylabel('y');
title('Van der Pol limit cycles and Fréchet mean (no phase alignment)');
legendStrings = [arrayfun(@(a) sprintf('a=%.2f',a), a_vals,'UniformOutput',false), {'Fréchet mean'}];
legend(legendStrings);

%% --------- Helper function: Arc-length resampling ---------
function curveOut = resampleArcLength(curveIn, M)
    diffs = diff(curveIn,1,1);
    segLen = sqrt(sum(diffs.^2,2));
    cumLen = [0; cumsum(segLen)];
    t_new = linspace(0, cumLen(end), M);
    x_new = interp1(cumLen, curveIn(:,1), t_new);
    y_new = interp1(cumLen, curveIn(:,2), t_new);
    curveOut = [x_new(:), y_new(:)];
end

end