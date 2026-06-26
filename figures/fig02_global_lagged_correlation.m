%% ================================================================
%  FIGURE 2: GLOBAL LAGGED CORRELATION ANALYSIS
%
%  PURPOSE:
%  Compute and visualise globally averaged (all-MLT) lagged
%  correlations between EEP and the two solar wind coupling
%  functions (Kan-Lee KL and epsilon) for each energy channel.
%
%  OUTPUT:
%  3-panel figure (one row per energy channel)
%    - KL and epsilon lagged correlation curves
%    - Peak-lag markers and annotations
%
%  ---------------------------------------------------------------
%  REQUIRED WORKSPACE VARIABLES (loaded via load_data.m)
%
%  Flux       : [Nt x NL x NMLT x NE]  electron flux 4D array
%  KLz        : standardised Kan-Lee coupling function  [Nt x 1]
%  EPSz       : standardised epsilon coupling function  [Nt x 1]
%  L_use      : L-shell index vector (L = 4–6)
%  lags       : lag vector (0:maxLag hours)
%  energyLabels : cell array of energy channel strings
%
%  Note: KLz, EPSz, L_use, lags, and energyLabels are defined in
%  fig01_MLT_lagged_correlation.m. Run that script first, or
%  redefine them in load_data.m.
%
% ================================================================

%% ================================================================
% GLOBAL MLT AVERAGE
% ================================================================

MLT_use = find(MLTvals >= 0 & MLTvals < 24);   % all MLT bins

nE = size(Flux, 4);

%% ================================================================
% PREALLOCATE
% ================================================================

R_KL  = nan(nE, length(lags));
R_EPS = nan(nE, length(lags));

%% ================================================================
% COMPUTE GLOBAL LAGGED CORRELATIONS
% ================================================================

for e = 1:nE

    % Global flux: average over all MLT and L = 4–6, then log10
    J = squeeze(mean( ...
        mean(Flux(:, L_use, MLT_use, e), 2, 'omitnan'), ...
        3, 'omitnan'));

    J = log10(J);
    J(~isfinite(J)) = NaN;

    for L = 1:length(lags)

        lag = lags(L);

        X1 = KLz(1:end-lag);
        X2 = EPSz(1:end-lag);
        Y  = J(1+lag:end);

        idx1 = isfinite(X1) & isfinite(Y);
        idx2 = isfinite(X2) & isfinite(Y);

        if sum(idx1) > 10
            R_KL(e,L) = corr(X1(idx1), Y(idx1));
        end

        if sum(idx2) > 10
            R_EPS(e,L) = corr(X2(idx2), Y(idx2));
        end

    end
end

%% ================================================================
% PLOTTING
% ================================================================

figure('Color', 'w', 'Position', [100 100 900 900]);

for e = 1:nE

    subplot(nE, 1, e)
    hold on

    % Correlation curves
    plot(lags, R_KL(e,:),  'LineWidth', 2);
    plot(lags, R_EPS(e,:), 'LineWidth', 2);

    % Peak markers
    [pk1, idx1] = max(R_KL(e,:));
    [pk2, idx2] = max(R_EPS(e,:));

    plot(lags(idx1), pk1, 'o', 'MarkerSize', 7, 'LineWidth', 1.5)
    plot(lags(idx2), pk2, 's', 'MarkerSize', 7, 'LineWidth', 1.5)

    ylim([0 0.7])

    % In-panel annotations
    text(lags(idx1)+0.5, pk1, ...
        sprintf('KL: r=%.2f @ %dh', pk1, lags(idx1)), 'FontSize', 10);
    text(lags(idx2)+0.5, pk2, ...
        sprintf('\\epsilon: r=%.2f @ %dh', pk2, lags(idx2)), 'FontSize', 10);

    xlabel('Lag [Hours]')
    ylabel('Correlation')

    if e == 1
        legend('KL', '\epsilon', 'Location', 'northeast')
    end

    if e == 2
        legend('', '', 'Peak KL', 'Peak \epsilon', 'Location', 'southeast');
    end

    grid on; box on;

end
