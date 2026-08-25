%% ================================================================
%  FIGURE 3 + TABLES 1, 2, 3: INTEGRATION WINDOW ANALYSIS
%
%  PRODUCES:
%   Figure 3  – Correlation vs. integration window (1–240 h)
%   Table 1   – Instantaneous Pearson correlations (KL and ε)
%   Table 2   – Peak integration window and peak r per channel
%   Table 3   – Multi-model regression summary + VIF
%
%  ---------------------------------------------------------------
%  REQUIRED WORKSPACE (from load_data.m)
%
%  fluxdata_2009_2019H_4D_new : [nVar x nTime x nMLT x nL]
%  KL_z, EPS_z                : standardised coupling functions
%  KL, EPS                    : raw coupling functions
%  logJp30, logJp100, logJp300: log10 precipitating flux series
%  Lidx                       : L = 4–6 index vector
%
% ================================================================

%% ================================================================
%% 1. EXTRACT GLOBAL FLUX SERIES (L = 4–6, all MLT)
%% ================================================================

mep0e1_3D = squeeze(fluxdata_2009_2019H_4D_new(7, :, :, :));
mep0e2_3D = squeeze(fluxdata_2009_2019H_4D_new(8, :, :, :));
mep0e3_3D = squeeze(fluxdata_2009_2019H_4D_new(9, :, :, :));

mep0e1_ts = squeeze(mean(mean(mep0e1_3D(:, :, Lidx), 2, 'omitmissing'), 3, 'omitmissing'));
mep0e2_ts = squeeze(mean(mean(mep0e2_3D(:, :, Lidx), 2, 'omitmissing'), 3, 'omitmissing'));
mep0e3_ts = squeeze(mean(mean(mep0e3_3D(:, :, Lidx), 2, 'omitmissing'), 3, 'omitmissing'));

% Quality screen: E1 below instrument threshold set to NaN
mep0e1_ts(mep0e1_ts <= 100) = NaN;

%% ================================================================
%% 2. LOG-TRANSFORM
%% ================================================================

logJp30  = log10(mep0e1_ts);
logJp100 = log10(mep0e2_ts);
logJp300 = log10(mep0e3_ts);

fluxSet = {logJp30, logJp100, logJp300};

%% ================================================================
%% TABLE 1: INSTANTANEOUS PEARSON CORRELATIONS
%% ================================================================

Table1 = table( ...
    [corr(logJp30,  KL_z, 'rows', 'complete'); ...
     corr(logJp100, KL_z, 'rows', 'complete'); ...
     corr(logJp300, KL_z, 'rows', 'complete')], ...
    [corr(logJp30,  EPS_z, 'rows', 'complete'); ...
     corr(logJp100, EPS_z, 'rows', 'complete'); ...
     corr(logJp300, EPS_z, 'rows', 'complete')], ...
    [corr(KL_z, EPS_z, 'rows', 'complete'); NaN; NaN], ...
    'VariableNames', {'r_KL', 'r_EPS', 'r_KL_EPS'}, ...
    'RowNames', {'>30 keV', '>100 keV', '>300 keV'});

disp('=== TABLE 1: Instantaneous correlations ===');
disp(Table1);

%% ================================================================
%% FIGURE 3 + TABLE 2: INTEGRATION WINDOW SWEEP
%% ================================================================

windows = 1:240;
corrKL  = zeros(3, numel(windows));
corrEPS = zeros(3, numel(windows));

for w = 1:numel(windows)
    T      = windows(w);
    KLint  = movsum(KL_z,  [T 0], 'omitnan');
    EPSint = movsum(EPS_z, [T 0], 'omitnan');

    for ch = 1:3
        corrKL(ch,w)  = corr(fluxSet{ch}, KLint,  'Rows', 'complete');
        corrEPS(ch,w) = corr(fluxSet{ch}, EPSint, 'Rows', 'complete');
    end
end

% --- Figure 3 ---
figure('Color', 'w')

chTitles = {'>30 keV [el.cm^{-2}s^{-1}sr^{-1}]', ...
            '>100 keV [el.cm^{-2}s^{-1}sr^{-1}]', ...
            '>300 keV [el.cm^{-2}s^{-1}sr^{-1}]'};

for ch = 1:3
    subplot(3, 1, ch);
    plot(windows, corrKL(ch,:),  'LineWidth', 2); hold on;
    plot(windows, corrEPS(ch,:), 'LineWidth', 2);
    ylabel('Correlation r');
    title(chTitles{ch});
    if ch == 3
        ylim([0 0.60]);
        xlabel('Integration Window [Hour]');
    else
        ylim([0 0.72]);
        xticklabels('');
    end
    if ch == 1
        legend('Kan–Lee', '\epsilon', ...
            'Location', 'northeast', 'Orientation', 'horizontal');
    end
    grid on; box on;
end

% --- Table 2 ---
peakVals = zeros(3, 4);
for ch = 1:3
    [pk1, idx1] = max(corrKL(ch,:));
    [pk2, idx2] = max(corrEPS(ch,:));
    peakVals(ch,:) = [pk1 windows(idx1) pk2 windows(idx2)];
end

Table2 = array2table(peakVals, ...
    'VariableNames', {'Peak_r_KL', 'PeakWindow_KL', 'Peak_r_EPS', 'PeakWindow_EPS'}, ...
    'RowNames', {'>30 keV', '>100 keV', '>300 keV'});

disp('=== TABLE 2: Peak integration window correlations ===');
disp(Table2);

%% ================================================================
%% TABLE 3: MULTI-MODEL REGRESSION + VIF
%%
%%  Integration windows from Table 2 peak windows:
%%    lagKL  = 6 h (>30 keV optimal)
%%    lagEPS = 15 h (>30 keV optimal)
%%  Adjust these values to match Table 2 output.
%% ================================================================

lagKL  = 6;
lagEPS = 15;

KL_int  = movmean(KL, [lagKL  0], 'omitmissing');
EPS_int = movmean(EPS, [lagEPS 0], 'omitmissing');

% Standardise integrated predictors
KL_inth  = (KL_int  - nanmean(KL_int))  ./ nanstd(KL_int);
EPS_inth = (EPS_int - nanmean(EPS_int)) ./ nanstd(EPS_int);

% Standardise response
Y_z = (logJp30 - nanmean(logJp30)) ./ nanstd(logJp30);

Xreg = [KL_inth EPS_inth];

% Z-scored models — for Table 3 standardised coefficients
mdl1 = fitlm(KL_inth,  Y_z);
mdl2 = fitlm(EPS_inth, Y_z);
mdl3 = fitlm(Xreg,     Y_z);

% Log-transformed but unstandardised models — for Figures 4 & 5
% Predictors: log10(integrated CF + offset); Response: logJp30
offsetKL  = 1e-3;
offsetEPS = 1e-3;

KLlog  = log10(KL_int  + offsetKL);
EPSlog = log10(EPS_int + offsetEPS);
Y      = logJp30;

Xreg_log = [KLlog EPSlog];

valid_log = isfinite(KLlog) & isfinite(EPSlog) & isfinite(Y);

mdl1_log = fitlm(KLlog(valid_log),               Y(valid_log));
mdl2_log = fitlm(EPSlog(valid_log),              Y(valid_log));
mdl3_log = fitlm(Xreg_log(valid_log,:),          Y(valid_log));

% Coefficients
a_KL1  = mdl1.Coefficients.Estimate(2);
a_EPS2 = mdl2.Coefficients.Estimate(2);
a_KL3  = mdl3.Coefficients.Estimate(2);
a_EPS3 = mdl3.Coefficients.Estimate(3);

p_KL1  = mdl1.Coefficients.pValue(2);
p_EPS2 = mdl2.Coefficients.pValue(2);
p_KL3  = mdl3.Coefficients.pValue(2);
p_EPS3 = mdl3.Coefficients.pValue(3);

R2_kl     = mdl1.Rsquared.Ordinary;
R2_eps    = mdl2.Rsquared.Ordinary;
R2_kl_eps = mdl3.Rsquared.Ordinary;

% Inter-predictor correlation
r_KL_EPS = corr(Xreg, 'rows', 'pairwise');
r_KL_EPS = r_KL_EPS(1,2);

% VIF (combined model)
nPred = 2;
VIF   = zeros(nPred, 1);
for j = 1:nPred
    yj    = Xreg(:, j);
    Xj    = Xreg(:, setdiff(1:nPred, j));
    mdlj  = fitlm(Xj, yj);
    VIF(j) = 1 / (1 - mdlj.Rsquared.Ordinary);
end

Table3 = table( ...
    {'KanLee'; 'Epsilon'; 'Kan-Lee + Epsilon'}, ...
    [R2_kl;   R2_eps;   R2_kl_eps], ...
    [a_KL1;   NaN;      a_KL3], ...
    [NaN;     a_EPS2;   a_EPS3], ...
    [p_KL1;   NaN;      p_KL3], ...
    [NaN;     p_EPS2;   p_EPS3], ...
    [NaN;     NaN;      VIF(1)], ...
    [NaN;     NaN;      VIF(2)], ...
    'VariableNames', {'Predictor','R2','Beta_KL','Beta_EPS', ...
                      'p_KL','p_EPS','VIF_KL','VIF_EPS'});

fprintf('\nr(KL, epsilon) = %.3f\n', r_KL_EPS);
disp('=== TABLE 3: Regression coefficients + VIF ===');
disp(Table3);
