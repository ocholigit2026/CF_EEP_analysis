%% ================================================================
%  utils/selectPeakEvents.m
%
%  From candidate event indices, select one representative peak
%  per cluster, where a cluster is a group of indices separated
%  by fewer than minSep steps.
%
%  INPUT
%    candidates : [N x 1]  sorted integer indices
%    series     : [Nt x 1]  time series (used to score peaks)
%    minSep     : integer   minimum inter-event separation
%
%  OUTPUT
%    selected : [Npeaks x 1]  one index per cluster (the maximum)
%
% ================================================================

function selected = selectPeakEvents(candidates, series, minSep)

    if isempty(candidates)
        selected = [];
        return
    end

    clusterID = cumsum([1; diff(candidates) > minSep]);

    selected = [];

    for k = 1:max(clusterID)
        members = candidates(clusterID == k);
        [~, idxMax] = max(series(members));
        selected(end+1) = members(idxMax);  %#ok<AGROW>
    end

    selected = selected(:);

end
