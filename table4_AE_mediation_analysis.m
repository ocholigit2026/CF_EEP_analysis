%% ================================================================
%  TABLE 4: AE MEDIATION ANALYSIS
%
%  PURPOSE:
%  Quantify the relative roles of external solar wind forcing (KL)
%  and internal magnetospheric activity (AE index) in driving EEP
%  variability, under both instantaneous and integrated KL forcing.
%
%  MODELS FITTED (per energy channel):
%    A:  EEP ~ ΣKL   and   EEP ~ KL (instantaneous)
%    B:  EEP ~ AE
%    C:  EEP ~ ΣKL + AE   and   EEP ~ KL + AE
%
%  Results in parentheses in the manuscript correspond to
%  instantaneous KL values.
%
%  ---------------------------------------------------------------
%  REQUIRED WORKSPACE (from load_data.m)
%
%  cpl_swdata : timetable with fields Ekl, AE,
%               flux30_0deg, flux100_0deg, flux300_0deg
%
% ================================================================

%% ================================================================
%% SETTINGS
%% ================================================================

intWindow    = 6;    % integration window (hours) — from Table 2
energyLabels = {'>30 keV', '>100 keV', '>300 keV'};
fluxFields   = {'flux30_0deg', 'flux100_0deg', 'flux300_0deg'};
nEnergy      = 3;

summaryRows  = {};

%% ================================================================
%% MAIN LOOP OVER ENERGY CHANNELS
%% ================================================================

for e = 1:nEnergy

    %% ----------------------------------------------------------
    %% 1. EXTRACT AND CLEAN
    %% ----------------------------------------------------------

    x1 = cpl_swdata.Ekl(:);
    x2 = cpl_swdata.AE(:);
    y  = cpl_swdata.(fluxFields{e})(:);

    valid = isfinite(x1) & isfinite(x2) & isfinite(y);
    x1 = x1(valid);
    x2 = x2(valid);
    y  = log10(y(valid));

    %% ----------------------------------------------------------
    %% 2. STANDARDISE
    %% ----------------------------------------------------------

    zAE  = (x2 - nanmean(x2)) ./ nanstd(x2);
    zY   = (y  - nanmean(y))  ./ nanstd(y);

    KL_int  = movmean(x1, [intWindow 0], 'omitmissing');
    zKL_int = (KL_int - nanmean(KL_int)) ./ nanstd(KL_int);
    zKL_ins = (x1     - nanmean(x1))     ./ nanstd(x1);

    %% ----------------------------------------------------------
    %% 3. FIT MODELS — INTEGRATED FORCING
    %% ----------------------------------------------------------

    mdl_KLint    = fitlm(zKL_int,           zY);
    mdl_AE       = fitlm(zAE,               zY);
    mdl_BOTH_int = fitlm([zKL_int zAE],     zY, 'VarNames', {'KL','AE','EEP'});

    R2_KLint    = mdl_KLint.Rsquared.Ordinary;
    beta_KLint  = mdl_KLint.Coefficients.Estimate(2);

    R2_AE       = mdl_AE.Rsquared.Ordinary;
    beta_AE     = mdl_AE.Coefficients.Estimate(2);

    R2_BOTH_int      = mdl_BOTH_int.Rsquared.Ordinary;
    beta_KL_both_int = mdl_BOTH_int.Coefficients.Estimate(2);
    beta_AE_both_int = mdl_BOTH_int.Coefficients.Estimate(3);

    %% ----------------------------------------------------------
    %% 4. FIT MODELS — INSTANTANEOUS FORCING
    %% ----------------------------------------------------------

    mdl_KLins    = fitlm(zKL_ins,           zY);
    mdl_BOTH_ins = fitlm([zKL_ins zAE],     zY, 'VarNames', {'KL','AE','EEP'});

    R2_KLins    = mdl_KLins.Rsquared.Ordinary;
    beta_KLins  = mdl_KLins.Coefficients.Estimate(2);

    R2_BOTH_ins      = mdl_BOTH_ins.Rsquared.Ordinary;
    beta_KL_both_ins = mdl_BOTH_ins.Coefficients.Estimate(2);
    beta_AE_both_ins = mdl_BOTH_ins.Coefficients.Estimate(3);

    %% ----------------------------------------------------------
    %% 5. PARTIAL CORRELATIONS (integrated forcing)
    %% ----------------------------------------------------------

    [r_KL_partial, p_KL_partial] = partialcorr(zKL_int, zY, zAE);
    [r_AE_partial, p_AE_partial] = partialcorr(zAE,     zY, zKL_int);

    %% ----------------------------------------------------------
    %% 6. PRINT
    %% ----------------------------------------------------------

    fprintf('\n====================================\n');
    fprintf('%s\n', energyLabels{e});
    fprintf('====================================\n');

    fprintf('\n--- Integrated ΣKL (window = %d h) ---\n', intWindow);
    fprintf('ΣKL only:    R2 = %.3f  beta = %.3f\n', R2_KLint, beta_KLint);
    fprintf('AE only:     R2 = %.3f  beta = %.3f\n', R2_AE,    beta_AE);
    fprintf('ΣKL + AE:   R2 = %.3f  betaKL = %.3f  betaAE = %.3f\n', ...
        R2_BOTH_int, beta_KL_both_int, beta_AE_both_int);
    fprintf('DeltaR2 adding AE to ΣKL : %.3f\n', R2_BOTH_int - R2_KLint);
    fprintf('DeltaR2 adding ΣKL to AE : %.3f\n', R2_BOTH_int - R2_AE);

    fprintf('\n--- Instantaneous KL ---\n');
    fprintf('KL only:     R2 = %.3f  beta = %.3f\n', R2_KLins, beta_KLins);
    fprintf('KL + AE:    R2 = %.3f  betaKL = %.3f  betaAE = %.3f\n', ...
        R2_BOTH_ins, beta_KL_both_ins, beta_AE_both_ins);

    fprintf('\n--- Partial correlations (integrated forcing) ---\n');
    fprintf('KL | AE:  r = %.3f  p = %.3e\n', r_KL_partial, p_KL_partial);
    fprintf('AE | KL:  r = %.3f  p = %.3e\n', r_AE_partial, p_AE_partial);

    %% ----------------------------------------------------------
    %% 7. COLLECT
    %% ----------------------------------------------------------

    summaryRows{end+1} = { ...
        energyLabels{e}, ...
        R2_KLint,  beta_KLint, ...
        R2_KLins,  beta_KLins, ...
        R2_AE,     beta_AE, ...
        R2_BOTH_int, beta_KL_both_int, beta_AE_both_int, ...
        R2_BOTH_ins, beta_KL_both_ins, beta_AE_both_ins};  %#ok<AGROW>

end

%% ================================================================
%% TABLE 4 SUMMARY
%% ================================================================

fprintf('\n\n=== TABLE 4: AE Mediation Summary ===\n');
fprintf('%-12s  %6s %6s  %6s %6s  %6s %6s  %8s %8s %8s  %8s %8s %8s\n', ...
    'Channel', ...
    'R2_int','b_int','R2_ins','b_ins', ...
    'R2_AE','b_AE', ...
    'R2_b_int','bKL_int','bAE_int', ...
    'R2_b_ins','bKL_ins','bAE_ins');

for row = summaryRows
    r = row{1};
    fprintf('%-12s  %6.3f %6.3f  %6.3f %6.3f  %6.3f %6.3f  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f\n', ...
        r{1}, r{2}, r{3}, r{4}, r{5}, r{6}, r{7}, ...
        r{8}, r{9}, r{10}, r{11}, r{12}, r{13});
end
