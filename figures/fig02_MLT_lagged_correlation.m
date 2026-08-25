%% ================================================================
%  FIGURE 2: MLT-DEPENDENT CF-EEP LAGGED CORRELATIONS
%
%  4 x 3 panel figure (>30 keV only)
%    Rows    = MLT sectors (Night, Dawn, Day, Dusk)
%    Columns = L-shell ranges (L=2-4, L=4-6, L=6-10)
%
%  Particle data layout: [nVar x nTime x nMLT x nL]
%  Positive lag = solar wind coupling leads EEP.
%
%  Side output: coverageTable, MLT_sector_coverage_statistics.csv
%
%  ---------------------------------------------------------------
%  REQUIRED WORKSPACE (from load_data.m)
%
%  fluxdata_2009_2019H_4D_new : [nVar x nTime x nL x nMLT]
%  KL_z, EPS_z                : standardised coupling functions
%  L_values, MLT_values       : spatial grids
%
% ================================================================

%% ----------------------------------------------------------------
%  1. SETTINGS  (must match fig01 if run in same session)
%% ----------------------------------------------------------------

maxLag = 24;
lags   = 0:maxLag;

L_values   = 1:0.25:15;
MLT_values = 0:23;

L_ranges = [2 4; 4 6; 6 10];
L_labels = {'L = 2--4', 'L = 4--6', 'L = 6--10'};
nLranges = size(L_ranges, 1);

energyLabels = {'>30 keV', '>100 keV', '>300 keV'};
nEnergy = 3;

%% ----------------------------------------------------------------
%  2. MLT SECTOR DEFINITIONS
%% ----------------------------------------------------------------

MLT_sector = false(4, numel(MLT_values));

MLT_sector(1,:) = (MLT_values >= 21) | (MLT_values < 3);   % Night
MLT_sector(2,:) = (MLT_values >= 3)  & (MLT_values < 9);   % Dawn
MLT_sector(3,:) = (MLT_values >= 9)  & (MLT_values < 15);  % Day
MLT_sector(4,:) = (MLT_values >= 15) & (MLT_values < 21);  % Dusk

sectorLabels = {'Night (21--03)', 'Dawn (03--09)', 'Day (09--15)', 'Dusk (15--21)'};
nSector = numel(sectorLabels);

%% ----------------------------------------------------------------
%  3. EXTRACT ELECTRON CHANNELS
%     fluxdata layout: nVar x nTime x nL x nMLT
%% ----------------------------------------------------------------

fluxdata = fluxdata_2009_2019H_4D_new;

mep0e1 = squeeze(fluxdata(7,:,:,:));
mep0e2 = squeeze(fluxdata(8,:,:,:));
mep0e3 = squeeze(fluxdata(9,:,:,:));

%% ----------------------------------------------------------------
%  4. PREALLOCATE
%% ----------------------------------------------------------------

R_KL_MLT  = nan(4, nLranges, numel(lags));
R_EPS_MLT = nan(4, nLranges, numel(lags));

peakR_KL_MLT   = nan(4, nLranges);
peakR_EPS_MLT  = nan(4, nLranges);
peakLag_KL_MLT  = nan(4, nLranges);
peakLag_EPS_MLT = nan(4, nLranges);

nUnderlyingValid = nan(nEnergy, nSector, nLranges);
nHourlyValid     = nan(nEnergy, nSector, nLranges);
nHourlyTotal     = nan(nEnergy, nSector, nLranges);
pctHourlyValid   = nan(nEnergy, nSector, nLranges);

%% ----------------------------------------------------------------
%  5. ANALYSIS: COVERAGE STATS (all energies) + CORRELATIONS (>30 keV)
%% ----------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf('FIGURE 2: MLT-DEPENDENT >30 keV CORRELATIONS\n');
fprintf('============================================================\n');

for e = 1:nEnergy

    if e == 1,     flux = mep0e1;
    elseif e == 2, flux = mep0e2;
    else,          flux = mep0e3;
    end

    for s = 1:4

        MLT_use = MLT_sector(s,:);

        for lr = 1:nLranges

            L_use = L_values >= L_ranges(lr,1) & L_values <= L_ranges(lr,2);

            % flux layout: [nTime x nMLT x nL]
            thisflux = flux(:, MLT_use, L_use);

            nUnderlyingValid(e,s,lr) = sum(isfinite(thisflux), 'all');

            J = squeeze(mean(mean(thisflux, 2, 'omitnan'), 3, 'omitnan'));

            validHourly = isfinite(J) & J > 0;
            nHourlyValid(e,s,lr)   = sum(validHourly);
            nHourlyTotal(e,s,lr)   = numel(J);
            pctHourlyValid(e,s,lr) = 100 * nHourlyValid(e,s,lr) / nHourlyTotal(e,s,lr);

            J(J <= 0) = NaN;
            J = log10(J);

            if e == 1
                for k = 1:numel(lags)

                    lag = lags(k);

                    xKL  = KL_z(1:end-lag);
                    xEPS = EPS_z(1:end-lag);
                    y    = J(1+lag:end);

                    valid = isfinite(xKL) & isfinite(y);
                    if sum(valid) >= 3
                        R_KL_MLT(s,lr,k) = corr(xKL(valid), y(valid), 'rows', 'complete');
                    end

                    valid = isfinite(xEPS) & isfinite(y);
                    if sum(valid) >= 3
                        R_EPS_MLT(s,lr,k) = corr(xEPS(valid), y(valid), 'rows', 'complete');
                    end

                end

                [peakR_KL_MLT(s,lr),  idxKL]  = max(R_KL_MLT(s,lr,:),  [], 'omitnan');
                [peakR_EPS_MLT(s,lr), idxEPS] = max(R_EPS_MLT(s,lr,:), [], 'omitnan');

                peakLag_KL_MLT(s,lr)  = lags(idxKL);
                peakLag_EPS_MLT(s,lr) = lags(idxEPS);

                fprintf('%s | %s | KL = %.3f (%dh) | EPS = %.3f (%dh)\n', ...
                    sectorLabels{s}, L_labels{lr}, ...
                    peakR_KL_MLT(s,lr), peakLag_KL_MLT(s,lr), ...
                    peakR_EPS_MLT(s,lr), peakLag_EPS_MLT(s,lr));
            end

        end
    end
end

%% ----------------------------------------------------------------
%  6. COVERAGE TABLE
%% ----------------------------------------------------------------

nRows = nEnergy * nSector * nLranges;

Energy = strings(nRows,1);  LRange = strings(nRows,1);
Sector = strings(nRows,1);
ValidHourly = zeros(nRows,1);  TotalHourly = zeros(nRows,1);
PercentValid = zeros(nRows,1); UnderlyingValid = zeros(nRows,1);

k = 0;
for e = 1:nEnergy
    for il = 1:nLranges
        for s = 1:nSector
            k = k + 1;
            Energy(k)          = energyLabels{e};
            LRange(k)          = L_labels{il};
            Sector(k)          = sectorLabels{s};
            ValidHourly(k)     = nHourlyValid(e,s,il);
            TotalHourly(k)     = nHourlyTotal(e,s,il);
            PercentValid(k)    = pctHourlyValid(e,s,il);
            UnderlyingValid(k) = nUnderlyingValid(e,s,il);
        end
    end
end

coverageTable = table(Energy, LRange, Sector, ValidHourly, TotalHourly, PercentValid, UnderlyingValid);

fprintf('\n============================================================\n');
fprintf('COMPLETE COVERAGE TABLE\n');
fprintf('============================================================\n');
disp(coverageTable);

writetable(coverageTable, 'MLT_sector_coverage_statistics.csv');

%% ----------------------------------------------------------------
%  7. FIGURE 2
%% ----------------------------------------------------------------

figure('Color', 'w', 'Position', [70 30 1200 1100]);

for s = 1:4
    for lr = 1:nLranges

        subplot(4, 3, (s-1)*3 + lr);
        hold on;

        RKL  = squeeze(R_KL_MLT(s,lr,:));
        REPS = squeeze(R_EPS_MLT(s,lr,:));

        plot(lags, RKL,  'LineWidth', 2);
        plot(lags, REPS, 'LineWidth', 2);

        plot(peakLag_KL_MLT(s,lr),  peakR_KL_MLT(s,lr),  'o', 'MarkerSize', 7, 'LineWidth', 1.5);
        plot(peakLag_EPS_MLT(s,lr), peakR_EPS_MLT(s,lr), 's', 'MarkerSize', 7, 'LineWidth', 1.5);

        xlim([0 maxLag]);
        ylim([0 0.7]);
        xticks(0:4:maxLag);
        yticks(0:0.1:0.7);
        grid on; box on;

        if s == 1
            title(L_labels{lr}, 'FontWeight', 'bold');
        end

        if lr == 1
            ylabel({sectorLabels{s}; 'Correlation'});
        end

        if s == 4 && lr == 2
            xlabel('Lag [Hour]');
        end

        text(0.04*maxLag, 0.63, sprintf('KL: %.2f @ %dh', peakR_KL_MLT(s,lr), peakLag_KL_MLT(s,lr)), 'FontSize', 8);
        text(0.04*maxLag, 0.56, sprintf('\\epsilon: %.2f @ %dh', peakR_EPS_MLT(s,lr), peakLag_EPS_MLT(s,lr)), 'FontSize', 8);

        if s == 1 && lr == 1
            legend('KL', '\epsilon', 'Location', 'northeast');
        end

    end
end

%% ----------------------------------------------------------------
%  8. SUMMARY OUTPUT
%% ----------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf('MLT PEAK LAGS (>30 keV)\n');
fprintf('============================================================\n');

for s = 1:4
    fprintf('\n%s\n', sectorLabels{s});
    for lr = 1:nLranges
        fprintf('%s: KL = %.3f @ %dh | EPS = %.3f @ %dh\n', ...
            L_labels{lr}, ...
            peakR_KL_MLT(s,lr), peakLag_KL_MLT(s,lr), ...
            peakR_EPS_MLT(s,lr), peakLag_EPS_MLT(s,lr));
    end
end

%% ----------------------------------------------------------------
%  9. SAVE
%% ----------------------------------------------------------------

save('Revised_Figure2_results.mat', ...
    'R_KL_MLT', 'R_EPS_MLT', ...
    'peakR_KL_MLT', 'peakR_EPS_MLT', ...
    'peakLag_KL_MLT', 'peakLag_EPS_MLT', ...
    'lags', 'L_values', 'MLT_values', 'L_ranges', 'sectorLabels');
