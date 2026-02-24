function new_outDir = validate_spike_counts(raw, outDir)
%VALIDATE_SPIKE_COUNTS Error out if bursts have different spike counts.

u = unique(raw.spikeCounts); % returns set of distinct values
if numel(u) ~= 1 
    % if u is not a single value
    fprintf('ATTENTION: bursts have different spike counts: %s \n', mat2str(raw.spikeCounts));
    new_outDir = fullfile(outDir, 'inconsistent_case');
else
    fprintf('All bursts validated: spikes=%d\n', u(1));
    new_outDir = fullfile(outDir, 'consistent_case');
end

end
