%% ================================================================
%  utils/classifyStormNonstorm.m
%
%  Classify a list of events as storm or nonstorm based on Dst.
%
%  INPUT
%    events    : [N x 1]  integer indices into Dst
%    Dst       : [Nt x 1] Dst time series
%    threshold : scalar   Dst threshold (e.g. -30 nT)
%    win       : integer  half-window around event to check (hours)
%
%  OUTPUT
%    stormIdx    : indices of storm-associated events
%    nonstormIdx : indices of nonstorm events
%
% ================================================================

function [stormIdx, nonstormIdx] = classifyStormNonstorm(events, Dst, threshold, win)

    stormIdx    = [];
    nonstormIdx = [];

    for i = 1:length(events)
        idx    = events(i);
        wStart = max(1, idx - win);
        wEnd   = min(length(Dst), idx + win);

        if min(Dst(wStart:wEnd)) <= threshold
            stormIdx(end+1) = idx;
        else
            nonstormIdx(end+1) = idx;
        end
    end

end
