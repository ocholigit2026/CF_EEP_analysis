%% ================================================================
%  TABLE 4: AE MEDIATION / PARTIAL CORRELATION ANALYSIS
%
%  PURPOSE:
%  Quantify the relative roles of:
%    (1) External solar wind forcing  – Kan-Lee (KL)
%    (2) Internal magnetospheric activity – AE index
%  in driving energetic electron precipitation (EEP).
%
%  MODELS FITTED:
%    A: EEP ~ KL           (external forcing alone)
%    B: EEP ~ AE           (internal activity alone)
%    C: EEP ~ KL + AE      (combined model)
%
%  DIAGNOSTICS:
%    - Ordinary R² for each model
%    - Standardised regression coefficients (beta)
%    - p-values
%    - Partial correlations  (KL | AE,  AE | KL)
%    - Incremental ΔR² when adding each predictor
%
%  ---------------------------------------------------------------
%  REQUIRED WORKSPACE VARIABLES (loaded via load_data.m)
%
%  cpl_swdata : timetable with fields  Ekl, AE, flux30_0deg
%
%  NOTE: The analysis uses a 7-hour cumulative KL window
%  (zKL_int = movsum(zKL, [7 0])).  Change KL_inp to zKL_ins
%  (instantaneous) for the direct-forcing variant.
%
% ================================================================

%% ================================================================
%% 1. EXTRACT AND CLEAN VARIABLES
%% ================================================================

x1 = cpl_swdata.Ekl(:);           % external forcing proxy (KL)
x2 = cpl_swdata.AE(:);            % internal activity proxy (AE)
y  = cpl_swdata.flux30_0deg(:);   % >30 keV precipitating flux

valid = isfinite(x1) & isfinite(x2) & isfinite(y);
x1 = x1(valid);
x2 = x2(valid);
y  = y(valid);

% Log-transform flux
y = log10(y);

%% ================================================================
%% 2. STANDARDISE
%% ================================================================

zKL = zscore(x1);
zAE = zscore(x2);
zY  = zscore(y);

% --- Forcing options ---
zKL_int = movsum(zKL, [7 0], 'omitmissing');   % 7-hour cumulative
zKL_ins = zKL;                                 % instantaneous

KL_inp = zKL_int;   % ← switch to zKL_ins for instantaneous variant

%% ================================================================
%% 3. MODEL A: EEP ~ KL
%% ================================================================

mdl_KL   = fitlm(KL_inp, zY);
R2_KL    = mdl_KL.Rsquared.Ordinary;
beta_KL  = mdl_KL.Coefficients.Estimate(2);
p_KL     = mdl_KL.Coefficients.pValue(2);

fprintf('\n====================================\n');
fprintf('MODEL A: EEP ~ KL\n');
fprintf('====================================\n');
fprintf('R^2      = %.3f\n', R2_KL);
fprintf('Beta(KL) = %.3f\n', beta_KL);
fprintf('p-value  = %.3e\n', p_KL);

%% ================================================================
%% 4. MODEL B: EEP ~ AE
%% ================================================================

mdl_AE   = fitlm(zAE, zY);
R2_AE    = mdl_AE.Rsquared.Ordinary;
beta_AE  = mdl_AE.Coefficients.Estimate(2);
p_AE     = mdl_AE.Coefficients.pValue(2);

fprintf('\n====================================\n');
fprintf('MODEL B: EEP ~ AE\n');
fprintf('====================================\n');
fprintf('R^2      = %.3f\n', R2_AE);
fprintf('Beta(AE) = %.3f\n', beta_AE);
fprintf('p-value  = %.3e\n', p_AE);

%% ================================================================
%% 5. MODEL C: EEP ~ KL + AE
%% ================================================================

X = [KL_inp zAE];

mdl_BOTH = fitlm(X, zY, 'VarNames', {'KL', 'AE', 'EEP'});
R2_BOTH  = mdl_BOTH.Rsquared.Ordinary;

beta_KL_both = mdl_BOTH.Coefficients.Estimate(2);
beta_AE_both = mdl_BOTH.Coefficients.Estimate(3);
p_KL_both    = mdl_BOTH.Coefficients.pValue(2);
p_AE_both    = mdl_BOTH.Coefficients.pValue(3);

fprintf('\n====================================\n');
fprintf('MODEL C: EEP ~ KL + AE\n');
fprintf('====================================\n');
fprintf('R^2           = %.3f\n', R2_BOTH);
fprintf('\nKL contribution:\n');
fprintf('Beta = %.3f\n', beta_KL_both);
fprintf('p    = %.3e\n', p_KL_both);
fprintf('\nAE contribution:\n');
fprintf('Beta = %.3f\n', beta_AE_both);
fprintf('p    = %.3e\n', p_AE_both);

%% ================================================================
%% 6. PARTIAL CORRELATIONS
%% ================================================================

[r_KL_partial, p_KL_partial] = partialcorr(KL_inp, zY, zAE);
[r_AE_partial, p_AE_partial] = partialcorr(zAE,    zY, KL_inp);

fprintf('\n====================================\n');
fprintf('PARTIAL CORRELATIONS\n');
fprintf('====================================\n');
fprintf('\nKL vs EEP (controlling AE):\n');
fprintf('r = %.3f\n', r_KL_partial);
fprintf('p = %.3e\n', p_KL_partial);
fprintf('\nAE vs EEP (controlling KL):\n');
fprintf('r = %.3f\n', r_AE_partial);
fprintf('p = %.3e\n', p_AE_partial);

%% ================================================================
%% 7. VARIANCE EXPLAINED (TABLE 4 SUMMARY)
%% ================================================================

fprintf('\n====================================\n');
fprintf('VARIANCE EXPLAINED  (TABLE 4)\n');
fprintf('====================================\n');
fprintf('KL alone         R^2 = %.3f\n', R2_KL);
fprintf('AE alone         R^2 = %.3f\n', R2_AE);
fprintf('KL + AE combined R^2 = %.3f\n', R2_BOTH);
fprintf('\nDelta R^2 adding AE to KL : %.3f\n', R2_BOTH - R2_KL);
fprintf('Delta R^2 adding KL to AE : %.3f\n', R2_BOTH - R2_AE);
