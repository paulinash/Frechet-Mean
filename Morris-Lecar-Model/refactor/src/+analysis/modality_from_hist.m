function out = modality_from_hist(x, opts)
%MODALITY_FROM_HIST Simple unimodal/multimodal heuristic via smoothed histogram peaks.
x = x(:);
nb = opts.numBins;
sigma = opts.smoothSigma;

% histogram
edges = linspace(min(x), max(x), nb+1);
counts = histcounts(x, edges);
centers = 0.5*(edges(1:end-1)+edges(2:end));

% smooth counts with a small Gaussian kernel
if sigma > 0
    kHalf = max(1, ceil(3*sigma));
    t = -kHalf:kHalf;
    g = exp(-(t.^2)/(2*sigma^2));
    g = g / sum(g);
    countsSm = conv(counts, g, 'same');
else
    countsSm = counts;
end

% Find interior peaks
peaks = [];
for i = 2:numel(countsSm)-1
    if countsSm(i) > countsSm(i-1) && countsSm(i) > countsSm(i+1)
        peaks(end+1) = i;
    end
end

% Find edge peaks
if numel(countsSm) >= 2
    if countsSm(1) > countsSm(2)
        peaks = [1 peaks];
    end
    if countsSm(end) > countsSm(end-1)
        peaks = [peaks numel(countsSm)];
    end
end


out = struct();
out.edges = edges;
out.centers = centers;
out.counts = counts;
out.countsSm = countsSm;
out.numPeaks = numel(peaks);
out.peakLocations = centers(peaks);

% Rule-of-thumb classification
if out.numPeaks <= 1
    out.class = "unimodal";
else
    out.class = "multimodal";
end

end
