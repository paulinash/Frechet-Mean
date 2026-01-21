%% i think this is unnecessary

function [dvec, pairs] = pairwise_distances(curves, maxN)
%PAIRWISE_DISTANCES Pairwise distances (subsample curves if N is large).
N = numel(curves);

if N > maxN
    idx = randperm(N, maxN);
    curvesSub = curves(idx);
else
    idx = 1:N;
    curvesSub = curves;
end

Ns = numel(curvesSub);
np = Ns*(Ns-1)/2;
dvec = zeros(np,1);
pairs = zeros(np,2);

k = 0;
for i = 1:Ns
    for j = i+1:Ns
        k = k + 1;
        dvec(k) = analysis.curve_distance(curvesSub{i}, curvesSub{j});
        pairs(k,:) = [idx(i), idx(j)];
    end
end

end
