function n = count_spikes(Vseg)
%COUNT_SPIKES Count prominent voltage peaks within a burst segment.
%
% Mirrors your original heuristic:
% - findpeaks with MinPeakProminence = 0.2
% - keep peaks with prominence >= max(0.25*max(prominence), 0.2)

[~, ~, ~, proms] = findpeaks(Vseg, 'MinPeakProminence', 0.2);

if isempty(proms)
    n = 0;
    return;
end

% consider only prominent spikes
prom_thr = max(0.25*max(proms), 0.2);
n = sum(proms >= prom_thr);
end
