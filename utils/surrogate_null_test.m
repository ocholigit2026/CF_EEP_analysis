%% ================================================================
%  SURROGATE NULL TEST: CIRCULAR-SHIFT ANALYSIS
%
%  PURPOSE:
%  Test whether the peak correlation obtained by optimizing the
%  temporal integration window in Table 2 exceeds the correlation
%  expected from chance temporal alignment.
%
%  The surrogate statistic is an exact replica of the Table 2
%  integration-window statistic:
%    - integration windows = 1:240 h
%    - zero lag only
%    - backward-looking movsum(CFz,[T 0],'omitnan')
%    - standardized KL and epsilon coupling functions
%    - log-transformed EEP
%    - Pearson correlation with complete observations
%    - peak statistic = maximum correlation across all 240 windows
%
%  For each surrogate realization, the coupling functions are
%  circularly shifted while the EEP series is left unchanged.
%  KL and epsilon are shifted by the same amount so their mutual
%  temporal relationship is preserved.
%
%  REQUIRED WORKSPACE:
%  Run load_data.m once before running this script.
%
%  Required variables:
%    fluxdata_2009_2019H_4D_new : [nVar x nTime x nMLT x nL]
%    cpl_swdata                 : timetable containing Ekl and eps
%
%  OUTPUTS:
%    windows    : integration windows used in Table 2
%    nSur       : number of surrogate realizations
%    R_obs      : observed correlations [3 x 2 x nW]
%    obsPeakR   : observed peak correlations [3 x 2]
%    obsPeakW   : corresponding optimal windows [3 x 2]
%    surPeakR   : surrogate peak correlations [nSur x 3 x 2]
%    nullMean   : mean surrogate peak correlation [3 x 2]
%    global95   : 95th percentile of surrogate peak correlations [3 x 2]
%    p_global   : global empirical p-values [3 x 2]
%
% ================================================================

%% ================================================================
%% 1. REQUIRED WORKSPACE
% ================================================================

requiredVars = { ...
    'fluxdata_2009_2019H_4D_new', ...
    'cpl_swdata'};

for k = 1:numel(requiredVars)

    if ~exist(requiredVars{k}, 'var')

        error(['Required workspace variable "%s" is missing. ' ...
               'Run load_data.m first.'], ...
               requiredVars{k});

    end

end

%% ================================================================
%% 2. SETTINGS
% ================================================================

% These settings exactly match the Table 2 integration analysis.
windows = 1:240;

% Number of circular-shift surrogate realizations.
nSur = 500;

% Minimum circular shift used for the surrogate construction.
%
% This is intentionally larger than the longest integration window.
% It prevents very small shifts from retaining substantial local
% temporal overlap between the original and shifted coupling series.
%
% This restriction affects only the surrogate construction; it does
% not affect the observed Table 2 statistic.
minShift = max(windows) + 1;

% Reproducible randomization.
rng(20260818, 'twister');

%% ================================================================
%% 3. EXTRACT GLOBAL EEP SERIES: L = 4–6, ALL MLT
% ================================================================

% L-shell grid used by the prepared MEPED product.
L_values = 1:0.25:15;

% Select L = 4–6.
Lidx = find(L_values >= 4 & L_values <= 6);

% Reference to the prepared 4-D satellite product:
%
%   variable x time x MLT x L
%
Flux = fluxdata_2009_2019H_4D_new;

% ---------------------------------------------------------------
% MEPED integral electron channels
% ---------------------------------------------------------------
%
% Variables 7, 8, and 9 correspond to:
%
%   variable 7 : >30 keV, 0-degree telescope
%   variable 8 : >100 keV, 0-degree telescope
%   variable 9 : >300 keV, 0-degree telescope
%
% The averaging below follows the same nested mean structure used
% in the Table 2 integration analysis:
%
%   first average over MLT,
%   then average over the selected L-shell bins.
%
mep0e1_3D = squeeze(Flux(7, :, :, :));
mep0e2_3D = squeeze(Flux(8, :, :, :));
mep0e3_3D = squeeze(Flux(9, :, :, :));

% ---------------------------------------------------------------
% Global hourly EEP series for L = 4–6
% ---------------------------------------------------------------

Jp30 = squeeze( ...
    mean( ...
        mean(mep0e1_3D(:, :, Lidx), 2, 'omitmissing'), ...
        3, 'omitmissing'));

Jp100 = squeeze( ...
    mean( ...
        mean(mep0e2_3D(:, :, Lidx), 2, 'omitmissing'), ...
        3, 'omitmissing'));

Jp300 = squeeze( ...
    mean( ...
        mean(mep0e3_3D(:, :, Lidx), 2, 'omitmissing'), ...
        3, 'omitmissing'));

% Force column-vector orientation.
Jp30  = Jp30(:);
Jp100 = Jp100(:);
Jp300 = Jp300(:);

%% ================================================================
%% 4. COUPLING FUNCTIONS AND STANDARDIZATION
% ================================================================

% Extract coupling functions from the synchronized timetable.
KL  = cpl_swdata.Ekl(:);
EPS = cpl_swdata.eps(:);

% Standardize exactly as in Table 2.
KLz = ( ...
    KL - mean(KL, 'omitnan')) ./ ...
    std(KL, 'omitnan');

EPSz = ( ...
    EPS - mean(EPS, 'omitnan')) ./ ...
    std(EPS, 'omitnan');

%% ================================================================
%% 5. LENGTH CHECK
% ================================================================

N = numel(KLz);

assert( ...
    numel(EPSz) == N, ...
    'KL and epsilon length mismatch.');

assert( ...
    numel(Jp30) == N, ...
    'KL and >30 keV EEP length mismatch.');

assert( ...
    numel(Jp100) == N, ...
    'KL and >100 keV EEP length mismatch.');

assert( ...
    numel(Jp300) == N, ...
    'KL and >300 keV EEP length mismatch.');

fprintf('\n');
fprintf('============================================================\n');
fprintf('CIRCULAR-SHIFT SURROGATE NULL TEST\n');
fprintf('============================================================\n');

fprintf('Hourly observations : %d\n', N);
fprintf('Integration windows : %d–%d h\n', ...
    min(windows), max(windows));
fprintf('Lag                 : 0 h only\n');
fprintf('Surrogates          : %d\n', nSur);
fprintf('Minimum shift       : %d h\n', minShift);

%% ================================================================
%% 6. PREPARE EEP EXACTLY AS IN TABLE 2
% ================================================================

% ---------------------------------------------------------------
% >30 keV quality screen
% ---------------------------------------------------------------
%
% The final Table 2 analysis applies the quality screen only to the
% >30 keV channel.
%
% Values <= 100 are treated as invalid.
%
Jp30(Jp30 <= 100) = NaN;

% ---------------------------------------------------------------
% Log transformation
% ---------------------------------------------------------------
%
% Apply the same log10 transformation used by Table 2.
%
logJp30  = log10(Jp30);
logJp100 = log10(Jp100);
logJp300 = log10(Jp300);

% Store the three EEP channels in a common cell array.
fluxSet = { ...
    logJp30, ...
    logJp100, ...
    logJp300};

energyLabels = { ...
    '>30 keV', ...
    '>100 keV', ...
    '>300 keV'};

CF_labels = { ...
    'KL', ...
    '\epsilon'};

%% ================================================================
%% 7. OBSERVED TABLE 2 STATISTIC
% ================================================================

% R_obs(e,c,w) contains the Pearson correlation for:
%
%   e = electron energy channel
%   c = coupling function
%   w = integration window
%
% The statistic is calculated at ZERO LAG only.
%
R_obs = nan(3, 2, numel(windows));

fprintf('\n');
fprintf('Calculating observed Table 2 statistic ...\n');

for iw = 1:numel(windows)

    T = windows(iw);

    % -----------------------------------------------------------
    % EXACT Table 2 integration convention
    % -----------------------------------------------------------
    %
    % The backward-looking window is implemented exactly as:
    %
    %   movsum(CFz,[T 0],'omitnan')
    %
    % Thus, for a given T, the current sample and T preceding
    % samples are included.
    %
    KLint = movsum( ...
        KLz, ...
        [T 0], ...
        'omitnan');

    EPSint = movsum( ...
        EPSz, ...
        [T 0], ...
        'omitnan');

    % -----------------------------------------------------------
    % Correlations for each energy channel
    % -----------------------------------------------------------

    for ch = 1:3

        R_obs(ch, 1, iw) = corr( ...
            fluxSet{ch}, ...
            KLint, ...
            'Rows', 'complete');

        R_obs(ch, 2, iw) = corr( ...
            fluxSet{ch}, ...
            EPSint, ...
            'Rows', 'complete');

    end

end

%% ================================================================
%% 8. OBSERVED PEAKS: EXACTLY AS TABLE 2
% ================================================================

% Peak correlation and corresponding integration window.
obsPeakR = nan(3, 2);
obsPeakW = nan(3, 2);

for ch = 1:3

    for c = 1:2

        rvec = squeeze(R_obs(ch, c, :));

        [obsPeakR(ch, c), idx] = max(rvec);

        obsPeakW(ch, c) = windows(idx);

    end

end

%% ================================================================
%% 9. VALID CIRCULAR SHIFTS
% ================================================================

% Only shifts greater than the maximum integration window are used.
%
% The shift is also kept away from the equivalent near-zero
% circular shift at the opposite end of the record.
%
possibleShifts = minShift:(N - minShift);

if isempty(possibleShifts)

    error([ ...
        'No valid circular shifts are available. ' ...
        'Reduce minShift or increase the record length.']);

end

%% ================================================================
%% 10. SURROGATE LOOP
% ================================================================

% Store only the maximum correlation obtained across the 240
% integration windows for each surrogate.
%
% Dimensions:
%
%   surrogate x energy channel x coupling function
%
%   [nSur x 3 x 2]
%
surPeakR = nan(nSur, 3, 2);

fprintf('\n');
fprintf('============================================================\n');
fprintf('RUNNING SURROGATE ANALYSIS\n');
fprintf('============================================================\n');

for s = 1:nSur

    % Progress message.
    if mod(s, 25) == 0 || s == 1

        fprintf( ...
            'Surrogate %d of %d\n', ...
            s, nSur);

    end

    % -----------------------------------------------------------
    % Select a random circular shift
    % -----------------------------------------------------------

    shiftAmt = possibleShifts( ...
        randi(numel(possibleShifts)));

    % -----------------------------------------------------------
    % Shift KL and epsilon by the SAME amount
    % -----------------------------------------------------------
    %
    % This preserves the internal temporal relationship between
    % the two coupling functions while destroying their alignment
    % with the EEP series.
    %
    KL_sur = circshift( ...
        KLz, ...
        shiftAmt);

    EPS_sur = circshift( ...
        EPSz, ...
        shiftAmt);

    % -----------------------------------------------------------
    % Initialize surrogate maxima
    % -----------------------------------------------------------

    peakKL  = -Inf(3, 1);
    peakEPS = -Inf(3, 1);

    % -----------------------------------------------------------
    % Reproduce the Table 2 window search
    % -----------------------------------------------------------

    for iw = 1:numel(windows)

        T = windows(iw);

        % EXACT Table 2 integration convention.
        KLint_sur = movsum( ...
            KL_sur, ...
            [T 0], ...
            'omitnan');

        EPSint_sur = movsum( ...
            EPS_sur, ...
            [T 0], ...
            'omitnan');

        % -------------------------------------------------------
        % Correlate with each EEP channel
        % -------------------------------------------------------

        for ch = 1:3

            rKL = corr( ...
                fluxSet{ch}, ...
                KLint_sur, ...
                'Rows', 'complete');

            rEPS = corr( ...
                fluxSet{ch}, ...
                EPSint_sur, ...
                'Rows', 'complete');

            % Update maximum KL correlation.
            if isfinite(rKL) && rKL > peakKL(ch)

                peakKL(ch) = rKL;

            end

            % Update maximum epsilon correlation.
            if isfinite(rEPS) && rEPS > peakEPS(ch)

                peakEPS(ch) = rEPS;

            end

        end

    end

    % -----------------------------------------------------------
    % Store surrogate maximum-over-window statistics
    % -----------------------------------------------------------

    surPeakR(s, :, 1) = peakKL;
    surPeakR(s, :, 2) = peakEPS;

end

fprintf('\n');
fprintf('Surrogate analysis complete.\n');

%% ================================================================
%% 11. NULL DISTRIBUTION STATISTICS
% ================================================================

% Mean of the surrogate maximum-over-window correlations.
nullMean = squeeze( ...
    mean(surPeakR, 1, 'omitnan'));

% 95th percentile of the surrogate maximum-over-window
% correlations.
global95 = squeeze( ...
    prctile(surPeakR, 95, 1));

%% ================================================================
%% 12. GLOBAL EMPIRICAL P-VALUES
% ================================================================

% The p-value tests the SAME maximum-over-window statistic
% reported in Table 2.
%
% The +1 correction is used for the finite number of surrogate
% realizations:
%
%       p = (1 + number of surrogate values >= observed value)
%           / (nSur + 1)
%
% With 500 surrogates, the smallest possible empirical p-value is
% approximately 0.002.
%
p_global = nan(3, 2);

for ch = 1:3

    for c = 1:2

        rsur = surPeakR(:, ch, c);

        % Retain only finite surrogate values.
        rsur = rsur(isfinite(rsur));

        if isempty(rsur)

            p_global(ch, c) = NaN;

        else

            p_global(ch, c) = ...
                (1 + sum( ...
                    rsur >= obsPeakR(ch, c))) ...
                / ...
                (numel(rsur) + 1);

        end

    end

end

%% ================================================================
%% 13. PRINT RESULTS
% ================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('GLOBAL NULL-TEST RESULTS\n');
fprintf('============================================================\n');

for ch = 1:3

    for c = 1:2

        fprintf('\n');
        fprintf('%s — %s\n', ...
            energyLabels{ch}, ...
            CF_labels{c});

        fprintf( ...
            'Observed Table 2 max R = %.4f\n', ...
            obsPeakR(ch, c));

        fprintf( ...
            'Optimal window         = %d h\n', ...
            obsPeakW(ch, c));

        fprintf( ...
            'Null mean max R        = %.4f\n', ...
            nullMean(ch, c));

        fprintf( ...
            'Null 95%% max R         = %.4f\n', ...
            global95(ch, c));

        fprintf( ...
            'Global p-value         = %.4f\n', ...
            p_global(ch, c));

    end

end

%% ================================================================
%% 14. SAVE RESULTS
% ================================================================

save( ...
    'EEP_integration_null_test.mat', ...
    'windows', ...
    'nSur', ...
    'R_obs', ...
    'surPeakR', ...
    'nullMean', ...
    'global95', ...
    'p_global', ...
    'obsPeakR', ...
    'obsPeakW');

fprintf('\n');
fprintf('============================================================\n');
fprintf('NULL TEST COMPLETE\n');
fprintf('Results saved: EEP_integration_null_test.mat\n');
fprintf('============================================================\n');
