%% ================================================================
%  CORRELATION AND REGRESSION MODULE
%
%  PRODUCES:
%   Figure 3  – Integration window length vs. correlation
%   Table  1  – Instantaneous Pearson correlations (KL and eps)
%   Table  2  – Peak integration window and peak r for each channel
%   Table  3  – Single-predictor regression + VIF (log-transformed)
%
%  NOTE ON TABLE 5 (reported in the manuscript):
%   Table 5 repeats the Table 3 analysis but uses log-transformed
%   versions of the instantaneous, lagged, AND integrated coupling
%   inputs to compute R-squared only.  Apply the same fitlm calls
%   below after replacing KLlog/EPSlog with your preferred
%   log-transformed predictor variant.
%
%  ---------------------------------------------------------------
%  REQUIRED WORKSPACE VARIABLES (loaded via load_data.m)
%
%  cpl_swdata : timetable with fields  Ekl, eps, flux*_0deg/90deg
%  t_hr       : datetime vector (hourly)
%
% ================================================================

%% ================================================================
%% 0. LOAD / EXTRACT PARTICLE DATA
%% ================================================================

% --- 0-degree telescope (precipitating) ---
mep0e1_3D = squeeze(fluxdata_2009_2019H_4D(7,  :, :, :));
mep0e2_3D = squeeze(fluxdata_2009_2019H_4D(8,  :, :, :));
mep0e3_3D = squeeze(fluxdata_2009_2019H_4D(9,  :, :, :));

% --- 90-degree telescope (trapped) ---
mep90e1_3D = squeeze(fluxdata_2009_2019H_4D(10, :, :, :));
mep90e2_3D = squeeze(fluxdata_2009_2019H_4D(11, :, :, :));
mep90e3_3D = squeeze(fluxdata_2009_2019H_4D(12, :, :, :));

% Helper: average over all MLT, then over L = 4–6 (indices 13:21)
avgFlux = @(x3D) ...
    mean(subsref(squeeze(mean(permute(x3D, [2 3 1]), 1, 'omitnan')), ...
            struct('type','()','subs',{{13:21, ':'}})), 1, 'omitnan')';

mep0e1_ts = avgFlux(mep0e1_3D);
mep0e2_ts = avgFlux(mep0e2_3D);
mep0e3_ts = avgFlux(mep0e3_3D);

mep90e1_ts = avgFlux(mep90e1_3D);
mep90e2_ts = avgFlux(mep90e2_3D);
mep90e3_ts = avgFlux(mep90e3_3D);

% Store in combined timetable
cpl_swdata.flux30_0deg  = mep0e1_ts;
cpl_swdata.flux100_0deg = mep0e2_ts;
cpl_swdata.flux300_0deg = mep0e3_ts;

cpl_swdata.flux30_90deg  = mep90e1_ts;
cpl_swdata.flux100_90deg = mep90e2_ts;
cpl_swdata.flux300_90deg = mep90e3_ts;

%% ================================================================
%% 1. PREPROCESSING
%% ================================================================

KL  = cpl_swdata.Ekl;
EPS = cpl_swdata.eps;

offset = 1e-3;   % small offset for log-transform stability

% --- Precipitating flux (0°) – quality screen for E1 ---
Jp30  = cpl_swdata.flux30_0deg;  Jp30(Jp30 <= 100) = NaN;
Jp100 = cpl_swdata.flux100_0deg;
Jp300 = cpl_swdata.flux300_0deg;

% --- Log-transform fluxes ---
logJp30  = log10(Jp30);
logJp100 = log10(Jp100);
logJp300 = log10(Jp300);

fluxSet = {logJp30, logJp100, logJp300};

% --- Standardise coupling proxies ---
KLz  = (KL(:)  - mean(KL,  'omitnan')) ./ std(KL,  'omitnan');
EPSz = (EPS(:) - mean(EPS, 'omitnan')) ./ std(EPS, 'omitnan');

%% ================================================================
%% TABLE 1: INSTANTANEOUS PEARSON CORRELATIONS
%% ================================================================

Table1 = table( ...
    [corr(logJp30,  KLz,  'rows', 'complete'); ...
     corr(logJp100, KLz,  'rows', 'complete'); ...
     corr(logJp300, KLz,  'rows', 'complete')], ...
    [corr(logJp30,  EPSz, 'rows', 'complete'); ...
     corr(logJp100, EPSz, 'rows', 'complete'); ...
     corr(logJp300, EPSz, 'rows', 'complete')], ...
    'VariableNames', {'r_KL', 'r_EPS'}, ...
    'RowNames',      {'>30 keV', '>100 keV', '>300 keV'});

disp('=== TABLE 1: Instantaneous correlations ===');
disp(Table1);

%% ================================================================
%% FIGURE 3 + TABLE 2: INTEGRATION WINDOW ANALYSIS
%% ================================================================

windows = 1:240;                     % window lengths 1–240 h
corrKL  = zeros(3, length(windows));
corrEPS = zeros(3, length(windows));

for w = 1:length(windows)

    T = windows(w);

    KLint  = movsum(KLz,  [T 0], 'omitnan');
    EPSint = movsum(EPSz, [T 0], 'omitnan');

    for ch = 1:3
        corrKL(ch,w)  = corr(fluxSet{ch}, KLint,  'Rows', 'complete');
        corrEPS(ch,w) = corr(fluxSet{ch}, EPSint, 'Rows', 'complete');
    end

end

% --- Figure 3 ---
figure('Color', 'w')

for ch = 1:3

    subplot(3, 1, ch);
    plot(windows, corrKL(ch,:),  'LineWidth', 2); hold on;
    plot(windows, corrEPS(ch,:), 'LineWidth', 2);
    ylabel('Correlation r');

    switch ch
        case 1
            title('>30 keV [el.cm^{-2}s^{-1}sr^{-1}]')
            ylim([0 0.72])
            legend('Kan–Lee', '\epsilon', ...
                'Location', 'northeast', 'Orientation', 'horizontal');
            xticklabels('')
        case 2
            title('>100 keV [el.cm^{-2}s^{-1}sr^{-1}]')
            ylim([0 0.72])
            xticklabels('')
        case 3
            title('>300 keV [el.cm^{-2}s^{-1}sr^{-1}]')
            ylim([0 0.60])
            xlabel('Integration Window [Hour]')
    end

    grid on; box on;

end

% --- Table 2: peak integration window ---
peakVals = zeros(3, 4);

for ch = 1:3
    [pk1, idx1] = max(corrKL(ch,:));
    [pk2, idx2] = max(corrEPS(ch,:));
    peakVals(ch,:) = [pk1 windows(idx1) pk2 windows(idx2)];
end

Table2 = array2table(peakVals, ...
    'VariableNames', {'Peak_r_KL', 'PeakWindow_KL', 'Peak_r_EPS', 'PeakWindow_EPS'}, ...
    'RowNames',      {'>30 keV', '>100 keV', '>300 keV'});

disp('=== TABLE 2: Integration window peaks ===');
disp(Table2);

%% ================================================================
%% TABLE 3: SINGLE-PREDICTOR REGRESSION + VIF
%%          (log-transformed, integrated forcing)
%%
%%  Uses peak integration windows from Table 2 (KL: 6 h, EPS: 15 h)
%%  as representative windows.  Adjust lagKL / lagEPS as needed.
%% ================================================================

lagKL  = 6;    % hours — adjust to peak window from Table 2
lagEPS = 15;

KL_int  = movmean(KL,  [lagKL  0], 'omitmissing');
EPS_int = movmean(EPS, [lagEPS 0], 'omitmissing');

% Log-transform integrated predictors
KLlog  = log10(KL_int  + offset);
EPSlog = log10(EPS_int + offset);

Y    = logJp30;
Xreg = [KLlog EPSlog];

% Individual models
mdl1 = fitlm(KLlog,  Y);
mdl2 = fitlm(EPSlog, Y);
mdl3 = fitlm(Xreg,   Y);

% Coefficients
a_KL1  = mdl1.Coefficients.Estimate(2);
a_EPS2 = mdl2.Coefficients.Estimate(2);

p_KL1  = mdl1.Coefficients.pValue(2);
p_EPS2 = mdl2.Coefficients.pValue(2);

R2_kl  = mdl1.Rsquared.Ordinary;
R2_eps = mdl2.Rsquared.Ordinary;

% VIF (variance inflation factor)
nPred = size(Xreg, 2);
VIF   = zeros(nPred, 1);

for j = 1:nPred
    yj     = Xreg(:, j);
    Xj     = Xreg(:, setdiff(1:nPred, j));
    mdl_j  = fitlm(Xj, yj);
    R2j    = mdl_j.Rsquared.Ordinary;
    VIF(j) = 1 / (1 - R2j);
end

Table3 = table( ...
    {'KanLee'; 'Epsilon'}, ...
    [R2_kl;  R2_eps], ...
    [a_KL1;  a_EPS2], ...
    [p_KL1;  p_EPS2], ...
    VIF, ...
    'VariableNames', {'Predictor', 'R2', 'Beta', 'pValue', 'VIF'});

disp('=== TABLE 3: Regression coefficients + VIF ===');
disp(Table3);
