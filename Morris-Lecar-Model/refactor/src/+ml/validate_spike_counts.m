function validate_spike_counts(raw)
%VALIDATE_SPIKE_COUNTS Error out if bursts have different spike counts.

u = unique(raw.spikeCounts); % returns set of distinct values
if numel(u) ~= 1 
    % if u is not a single value
    error('ABORT: bursts have different spike counts: %s', mat2str(raw.spikeCounts));
end
fprintf('All bursts validated: spikes=%d\n', u(1));
fprintf('Mean period = %.2f\n', mean(raw.periods));
end
