%% ================================================================
%  FIGURE 1: MLT-DEPENDENT LAGGED CORRELATION ANALYSIS
%
%  PURPOSE:
%  Examine whether spatial averaging obscures the lagged response
%  between EEP and solar wind coupling functions.
%
%  OUTPUT:
%  4 x 3 panel figure
%
%    ROWS    = MLT sectors (Night, Dawn, Day, Dusk)
%    COLUMNS = Electron energy channels (>30, >100, >300 keV)
%
%  Each panel contains:
%     - KL (Kan-Lee) lagged correlation curve
%     - epsilon lagged correlation curve
%     - peak lag/correlation markers and annotations
%
%  ---------------------------------------------------------------
%  REQUIRED WORKSPACE VARIABLES (loaded via load_data.m)
%
%  Flux          : [Nt x NL x NMLT x NE]  electron flux 4D array
%  cpl_swdata    : timetable with Ekl, eps fields
%
%  ---------------------------------------------------------------
%  DIMENSIONS
%
%  Nt   = number of hourly time steps
%  NL   = L-shell bins  (Lvals = 1:0.25:10)
%  NMLT = MLT bins
%  NE   = energy channels (3: >30, >100, >300 keV)
%
%  ---------------------------------------------------------------
%  AUTHOR NOTE
%
%  This implementation:
%    - averages flux over L = 4–6 and within broad MLT sectors
%    - log10-transforms the averaged flux
%    - computes Pearson lagged correlations (lag 0–24 h)
%    - identifies peak lag and correlation for each panel
%
% ================================================================

%% ================================================================
% USER SETTINGS
% ================================================================

maxLag = 24;          % maximum lag (hours)
lags   = 0:maxLag;   % positive lag: coupling leads precipitation

energyLabels = {'>30 keV', '>100 keV', '>300 keV'};

sectorLabels = { ...
    'Night (21-03)', ...
    'Dawn (03-09)', ...
    'Day (09-15)', ...
    'Dusk (15-21)'};

%% ================================================================
% DEFINE MLT SECTORS
% ================================================================

NMLT = size(Flux, 3);

Lvals   = 1:0.25:10;
L_use   = find(Lvals >= 4 & Lvals <= 6);

MLTvals = linspace(0, 24, NMLT+1);
MLTvals(end) = [];

sectorInds    = cell(4, 1);
sectorInds{1} = find(MLTvals >= 21 | MLTvals < 3);   % Night
sectorInds{2} = find(MLTvals >= 3  & MLTvals < 9);   % Dawn
sectorInds{3} = find(MLTvals >= 9  & MLTvals < 15);  % Day
sectorInds{4} = find(MLTvals >= 15 & MLTvals < 21);  % Dusk

%% ================================================================
% STANDARDIZE COUPLING FUNCTIONS
% ================================================================

KL  = cpl_swdata.Ekl(:);
EPS = cpl_swdata.eps(:);

KLz  = (KL  - mean(KL,  'omitnan')) ./ std(KL,  'omitmissing');
EPSz = (EPS - mean(EPS, 'omitnan')) ./ std(EPS, 'omitmissing');

%% ================================================================
% PREALLOCATE
% ================================================================

nSector = 4;
nEnergy = 3;

corrKL  = nan(nSector, nEnergy, length(lags));
corrEPS = nan(nSector, nEnergy, length(lags));

peakLag_KL  = nan(nSector, nEnergy);
peakLag_EPS = nan(nSector, nEnergy);
peakR_KL    = nan(nSector, nEnergy);
peakR_EPS   = nan(nSector, nEnergy);

%% ================================================================
% MAIN ANALYSIS: LAGGED CORRELATIONS
% ================================================================

for s = 1:nSector

    mltIdx = sectorInds{s};

    for e = 1:nEnergy

        % Sector- and L-averaged flux
        J = squeeze(mean(Flux(:, L_use, mltIdx, e), [2 3], 'omitnan'));

        % Log10-transform
        J = log10(J);

        for k = 1:length(lags)

            lag = lags(k);

            xKL  = KLz(1:end-lag);
            xEPS = EPSz(1:end-lag);
            y    = J(1+lag:end);

            valid1 = isfinite(xKL)  & isfinite(y);
            valid2 = isfinite(xEPS) & isfinite(y);

            if sum(valid1) > 10
                corrKL(s,e,k) = corr(xKL(valid1), y(valid1), 'rows', 'complete');
            end

            if sum(valid2) > 10
                corrEPS(s,e,k) = corr(xEPS(valid2), y(valid2), 'rows', 'complete');
            end

        end

        % Peak lag / correlation
        [peakR_KL(s,e),  idx1] = max(corrKL(s,e,:));
        [peakR_EPS(s,e), idx2] = max(corrEPS(s,e,:));

        peakLag_KL(s,e)  = lags(idx1);
        peakLag_EPS(s,e) = lags(idx2);

    end
end

%% ================================================================
% PLOTTING
% ================================================================

figure('Color', 'w', 'Position', [50 50 1500 1100]);

for s = 1:nSector
    for e = 1:nEnergy

        subplot(nSector, nEnergy, (s-1)*nEnergy + e)
        hold on

        % Correlation curves
        plot(lags, squeeze(corrKL(s,e,:)),  'LineWidth', 2);
        plot(lags, squeeze(corrEPS(s,e,:)), 'LineWidth', 2);

        % Peak markers
        plot(peakLag_KL(s,e),  peakR_KL(s,e),  'o', 'MarkerSize', 7, 'LineWidth', 1.5);
        plot(peakLag_EPS(s,e), peakR_EPS(s,e),  's', 'MarkerSize', 7, 'LineWidth', 1.5);

        % Formatting
        grid on; box on;

        if s == 4, xlabel('Lag [hours]'); end
        if e == 1, ylabel('Correlation');  end

        title(sprintf('%s | %s', sectorLabels{s}, energyLabels{e}))

        yline(0, 'k--')
        ylim([0 0.7])

        % Legend (shown only once per type)
        if s == 1 && e == 1
            legend('KL', '\epsilon', 'Location', 'northeast');
        end

        if s == 2 && e == 1
            legend('', '', 'Peak KL', 'Peak \epsilon', 'Location', 'southeast');
        end

        % In-panel annotations
        txt1 = sprintf('KL: r=%.2f @ %dh', peakR_KL(s,e), peakLag_KL(s,e));
        txt2 = sprintf('\\epsilon: r=%.2f @ %dh', peakR_EPS(s,e), peakLag_EPS(s,e));

        yl = ylim;
        text(0.02*maxLag, yl(2)-0.08*(yl(2)-yl(1)), txt1, 'FontSize', 8);
        text(0.02*maxLag, yl(2)-0.16*(yl(2)-yl(1)), txt2, 'FontSize', 8);

    end
end

%% ================================================================
% CONSOLE OUTPUT: PEAK CORRELATION TABLES
% ================================================================

disp('===================================================')
disp('PEAK KL CORRELATIONS')
disp('===================================================')
for s = 1:nSector
    fprintf('\n%s\n', sectorLabels{s});
    for e = 1:nEnergy
        fprintf('%s: r = %.3f at lag = %d h\n', energyLabels{e}, peakR_KL(s,e), peakLag_KL(s,e));
    end
end

disp(' ')
disp('===================================================')
disp('PEAK EPSILON CORRELATIONS')
disp('===================================================')
for s = 1:nSector
    fprintf('\n%s\n', sectorLabels{s});
    for e = 1:nEnergy
        fprintf('%s: r = %.3f at lag = %d h\n', energyLabels{e}, peakR_EPS(s,e), peakLag_EPS(s,e));
    end
end
