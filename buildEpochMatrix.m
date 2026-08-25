%% ================================================================
%  utils/buildEpochMatrix.m
%
%  Build a superposed-epoch matrix for a list of events.
%
%  INPUT
%    events : [Nevents x 1]  integer indices into series
%    series : [Nt x 1]       time series
%    win    : integer         half-window in steps (epoch = -win:win)
%
%  OUTPUT
%    M : [(2*win+1) x Nevents]
%        Columns outside the series boundaries are NaN.
%
% ================================================================

function M = buildEpochMatrix(events, series, win)

    tau = -win:win;
    nt  = length(tau);

    M = nan(nt, length(events));

    for i = 1:length(events)
        idx = events(i);

        if idx - win < 1 || idx + win > length(series)
            continue
        end

        M(:, i) = series(idx-win : idx+win);
    end

end
