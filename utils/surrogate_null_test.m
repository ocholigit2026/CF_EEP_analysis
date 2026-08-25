%% ================================================================
%  SURROGATE NULL TEST: CIRCULAR-SHIFT ANALYSIS
%
%  PURPOSE:
%  Test whether the correlation enhancement associated with temporal
%  integration exceeds what is expected from the autocorrelation
%  structure of each series alone, without actual CF–EEP alignment.
%
%  NULL HYPOTHESIS:
%  The observed CF–EEP relationship is explained by autocorrelation
%  and temporal smoothing, not by true temporal alignment.
%
%  METHOD:
%  Both KL and EPS are shifted by the same random circular amount,
%  preserving their mutual relationship. Surrogate correlations
%  are computed over all integration windows and lags.
%
%  OUTPUTS:
%    R_obs      [3 x 2 x nW x nLag]  observed correlations
%    R_sur      [nSur x 3 x 2 x nW x nLag]  surrogate correlations
%    surPeakR   [nSur x 3 x 2]  max surrogate R per realisation
%    obsPeakR   [3 x 2]  observed maximum R
%    p_global   [3 x 2]  global empirical p-values
%    global95   [3 x 2]  95th-percentile null thresholds
%
%  ---------------------------------------------------------------
%  REQUIRED WORKSPACE (from load_data.m)
%
%  fluxdata_2009_2019H_4D_new : [nVar x nTime x nMLT x nL]
%  KL_z, EPS_z                : standardised coupling functions
%  L_values, Lidx             : spatial grid and L=4-6 indices
%
% ================================================================

%% ================================================================
%% 1. SETTINGS
%% ================================================================

windows  = 1:49;
lags     = 0:24;
nSur     = 500;
minShift = max(windows) + max(lags) + 24;

rng(20260818, 'twister');

%% ================================================================
%% 2. EXTRACT GLOBAL EEP SERIES (L = 4–6, all MLT)
%% ================================================================

F30  = squeeze(mean(mean(fluxdata_2009_2019H_4D_new(7, :, :, Lidx), 3, 'omitnan'), 4, 'omitnan'));
F100 = squeeze(mean(mean(fluxdata_2009_2019H_4D_new(8, :, :, Lidx), 3, 'omitnan'), 4, 'omitnan'));
F300 = squeeze(mean(mean(fluxdata_2009_2019H_4D_new(9, :, :, Lidx), 3, 'omitnan'), 4, 'omitnan'));

F30  = F30(:);   F100 = F100(:);   F300 = F300(:);
F30(F30   <= 0) = NaN;
F100(F100 <= 0) = NaN;
F300(F300 <= 0) = NaN;

EEP = {log10(F30), log10(F100), log10(F300)};
energyLabels = {'>30 keV', '>100 keV', '>300 keV'};

CF_labels = {'KL', '\epsilon'};

nW   = numel(windows);
nLag = numel(lags);

fprintf('Hourly observations : %d\n', numel(F30));
fprintf('Min circular shift  : %d h\n', minShift);

%% ================================================================
%% 3. OBSERVED CORRELATIONS
%% ================================================================

R_obs = computeAllR({KL_z, EPS_z}, EEP, windows, lags);

%% ================================================================
%% 4. OBSERVED PEAK CORRELATIONS
%% ================================================================

obsPeakR   = nan(3, 2);
obsPeakW   = nan(3, 2);
obsPeakLag = nan(3, 2);

for e = 1:3
    for c = 1:2
        Rtmp = squeeze(R_obs(e, c, :, :));
        [obsPeakR(e,c), linIdx] = max(Rtmp(:), [], 'omitnan');
        [iw, ilag]              = ind2sub(size(Rtmp), linIdx);
        obsPeakW(e,c)           = windows(iw);
        obsPeakLag(e,c)         = lags(ilag);
    end
end

%% ================================================================
%% 5. SURROGATE LOOP
%% ================================================================

R_sur    = nan(nSur, 3, 2, nW, nLag);
surPeakR = nan(nSur, 3, 2);

possibleShifts = minShift:(numel(KL_z)-minShift);

fprintf('\nRunning %d surrogates ...\n', nSur);

for s = 1:nSur

    if mod(s, 50) == 0, fprintf('  Surrogate %d of %d\n', s, nSur); end

    shiftAmt = possibleShifts(randi(numel(possibleShifts)));
    CFs_sur  = {circshift(KL_z, shiftAmt), circshift(EPS_z, shiftAmt)};

    R_sur(s,:,:,:,:) = computeAllR(CFs_sur, EEP, windows, lags);

    for e = 1:3
        for c = 1:2
            Rtmp = squeeze(R_sur(s, e, c, :, :));
            surPeakR(s, e, c) = max(Rtmp(:), [], 'omitnan');
        end
    end

end

fprintf('Surrogate analysis complete.\n');

%% ================================================================
%% 6. NULL DISTRIBUTION AND P-VALUES
%% ================================================================

p_global = nan(3, 2);
global95 = nan(3, 2);

for e = 1:3
    for c = 1:2
        rsur = surPeakR(:, e, c);
        rsur = rsur(isfinite(rsur));
        p_global(e,c) = (1 + sum(rsur >= obsPeakR(e,c))) / (numel(rsur) + 1);
        global95(e,c) = prctile(rsur, 95);
    end
end

%% ================================================================
%% 7. PRINT RESULTS
%% ================================================================

fprintf('\n============================================================\n');
fprintf('GLOBAL NULL-TEST RESULTS\n');
fprintf('============================================================\n');

for e = 1:3
    for c = 1:2
        fprintf('\n%s — %s\n', energyLabels{e}, CF_labels{c});
        fprintf('Observed max R  = %.4f  (window=%dh, lag=%dh)\n', ...
            obsPeakR(e,c), obsPeakW(e,c), obsPeakLag(e,c));
        fprintf('Null 95%% max R  = %.4f\n', global95(e,c));
        fprintf('Global p-value  = %.4f\n', p_global(e,c));
    end
end

%% ================================================================
%% 8. FIGURES
%% ================================================================

% Figure 1: Null distributions vs observed peak R
figure('Color', 'w', 'Position', [100 100 1200 850]);
for e = 1:3
    subplot(3, 1, e); hold on;
    histogram(surPeakR(:,e,1), 25, 'Normalization', 'probability');
    histogram(surPeakR(:,e,2), 25, 'Normalization', 'probability');
    xline(obsPeakR(e,1), 'LineWidth', 2);
    xline(obsPeakR(e,2), 'LineWidth', 2);
    xline(global95(e,1), '--', 'LineWidth', 1.5);
    xline(global95(e,2), '--', 'LineWidth', 1.5);
    xlabel('Maximum correlation'); ylabel('Probability');
    title(energyLabels{e}); grid on; box on;
    if e == 1
        legend('KL null','\epsilon null','Observed KL','Observed \epsilon', ...
               'KL 95% null','\epsilon 95% null','Location','best');
    end
end
sgtitle('Observed peak correlations vs. circular-shift null distributions','FontWeight','bold');

% Figure 2: R(window, lag) heat maps
figure('Color', 'w', 'Position', [100 50 1300 950]);
for e = 1:3
    subplot(3,2,2*e-1);
    imagesc(windows, lags, squeeze(R_obs(e,1,:,:))');
    axis xy; xlabel('Integration window [h]'); ylabel('Lag [h]');
    title([energyLabels{e} ' — KL']);
    colorbar; caxis([0 0.7]); box on;

    subplot(3,2,2*e);
    imagesc(windows, lags, squeeze(R_obs(e,2,:,:))');
    axis xy; xlabel('Integration window [h]'); ylabel('Lag [h]');
    title([energyLabels{e} ' — \epsilon']);
    colorbar; caxis([0 0.7]); box on;
end
sgtitle('Observed correlation: integration window vs. lag','FontWeight','bold');

%% ================================================================
%% 9. SAVE
%% ================================================================

save('EEP_integration_null_test.mat', ...
    'windows', 'lags', 'nSur', ...
    'R_obs', 'R_sur', ...
    'p_global', 'obsPeakR', 'obsPeakW', 'obsPeakLag', ...
    'surPeakR', 'global95');

fprintf('\nResults saved: EEP_integration_null_test.mat\n');

%% ================================================================
%% LOCAL FUNCTION
%% ================================================================

function R = computeAllR(CFs, EEP, windows, lags)
    nW = numel(windows); nLag = numel(lags); nE = 3;
    R  = nan(nE, 2, nW, nLag);
    for e = 1:nE
        y = EEP{e};
        for c = 1:2
            CF = CFs{c};
            for iw = 1:nW
                CF_int = movsum(CF, [windows(iw)-1 0], 'omitnan');
                for ilag = 1:nLag
                    lag = lags(ilag);
                    if lag == 0
                        x = CF_int; yy = y;
                    else
                        x = CF_int(1:end-lag); yy = y(1+lag:end);
                    end
                    valid = isfinite(x) & isfinite(yy);
                    if sum(valid) >= 3
                        R(e,c,iw,ilag) = corr(x(valid), yy(valid), 'type', 'Pearson');
                    end
                end
            end
        end
    end
end
