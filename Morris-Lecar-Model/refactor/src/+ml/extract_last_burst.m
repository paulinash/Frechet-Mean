function [segT, segZ] = extract_last_burst(t, Z)
%EXTRACT_LAST_BURST Extract the last burst period using the slow variable y.
%
% - threshold based on y range
% - find downward crossings with negative gradient
% - take last two crossings as burst boundaries
%
% Inputs:
%   t: (T x 1) time after transient removal
%   Z: (T x 3) [V w y]

% Slow variable
y = Z(:,3);

y_min = min(y); y_max = max(y);
yth = y_min + 0.15*(y_max - y_min); % Low threshold on decline

dy = gradient(y, t); % Detect downward (!) crossings 

% Candidates: index where consecutive points (y(i),y(i+1)) cross threshold
downIdx = find((y(1:end-1) > yth) & (y(2:end) <= yth)) + 1;
% Choose index with actually decreasing gradient
downIdx = downIdx(downIdx <= numel(dy) & dy(downIdx) < -1e-4);

if numel(downIdx) < 2
    error('Not enough valid downward crossings to define a burst.');
end

i1 = downIdx(end-1);
i2 = downIdx(end);

segT = t(i1:i2);
segZ = Z(i1:i2,:);
end
