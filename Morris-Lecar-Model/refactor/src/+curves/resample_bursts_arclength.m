function [curvesOut, velsOut, info] = resample_bursts_arclength(raw, M)
%RESAMPLE_BURSTS_ARCLENGTH Resample burst segments to uniform normalized arc length.
%
% Input:
%   raw.segT{k} time, raw.segZ{k} states [V w y]
% Output:
%   curvesOut{k}  Mx3 curve sampled uniformly in normalized arc length
%   velsOut{k}    Mx3 velocities w.r.t. original time, resampled to same arc grid
%   info.sq       query grid in [0,1) (length M)

N = numel(raw.segT);
curvesOut = cell(N,1);
velsOut = cell(N,1);

sqFull = linspace(0,1,M+1);
sq = sqFull(1:end-1);
info = struct();
info.sqFull = sqFull;
info.sq = sq;

for k = 1:N
    t = raw.segT{k};
    Z = raw.segZ{k};

    % Velocities in time (finite difference via gradient)
    Vdot = gradient(Z(:,1), t);
    Wdot = gradient(Z(:,2), t);
    Ydot = gradient(Z(:,3), t);
    Vraw = [Vdot, Wdot, Ydot];

    % Arc length parameter s in [0,1]
    dZ = diff(Z,1,1);
    s = [0; cumsum(sqrt(sum(dZ.^2,2)))];
    sEnd = s(end);
    if sEnd <= 0
        error('Degenerate curve for k=%d: zero arc length.', k);
    end
    s = s / sEnd;

    % Interpolate to uniform arc grid
    curvesOut{k} = [ ...
        interp1(s, Z(:,1), sq, 'pchip')', ...
        interp1(s, Z(:,2), sq, 'pchip')', ...
        interp1(s, Z(:,3), sq, 'pchip')' ...
    ];
    velsOut{k} = interp1(s, Vraw, sq, 'pchip');
end
end
