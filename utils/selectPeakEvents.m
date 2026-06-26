%% ================================================================
%  utils/selectPeakEvents.m
%
%  From a list of candidate event indices (e.g. all time steps
%  exceeding a percentile threshold), select one representative
%  peak per cluster, where a cluster is defined as a group of
%  indices separated by fewer than minSep steps.
%
%  INPUT
%    candidates : [N x 1]  sorted integer indices of candidate events
%    series     : [Nt x 1]  the time series used to score peaks
%    minSep     : integer   minimum separation (time steps / hours)
%                           that defines cluster boundaries
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

    % Assign cluster IDs based on inter-event gaps
    clusterID = cumsum([1; diff(candidates) > minSep]);

    selected = [];

    for k = 1:max(clusterID)
        members = candidates(clusterID == k);

        % Keep the time step with the largest series value
        [~, idxMax] = max(series(members));
        selected(end+1) = members(idxMax);  %#ok<AGROW>
    end

    selected = selected(:);   % column vector

end
