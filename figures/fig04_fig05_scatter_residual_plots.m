%% ================================================================
%  FIGURES 4 & 5: PREDICTED vs OBSERVED SCATTER AND RESIDUAL PLOTS
%
%  FIGURE 4 - Predicted vs. observed scatter (3 panels):
%    - Density-coloured scatter
%    - 1:1 reference line
%    - Robust regression fit + 95% prediction interval band
%    - Binned mean +/- 1sigma error bars
%
%  FIGURE 5 - Residual diagnostics (3 panels):
%    - Scatter of residuals vs. predicted flux
%    - Binned mean +/- 1sigma (left axis)
%    - Percentage of observations per bin (right axis)
%
%  ---------------------------------------------------------------
%  REQUIRED WORKSPACE
%
%  mdl1_log, mdl2_log, mdl3_log : fitlm models (log-transformed,
%                                  unstandardised) from fig03_tables1_2_3_*.m
%  KLlog, EPSlog  : log10 of integrated coupling functions
%  Xreg_log       : [KLlog EPSlog]
%  valid_log      : logical index of finite observations
%  Y              : observed log10 precipitating flux (>30 keV)
%
%  NOTE: Run fig03_tables1_2_3_integration_window.m first.
%
% ================================================================

%% ----------------------------------------------------------------
%  1. GENERATE PREDICTIONS
%% ----------------------------------------------------------------

panelTitle = {'KL', '\epsilon', 'KL + \epsilon'};

% Preallocate prediction arrays.
pred_KL   = nan(size(Y));
pred_eps  = nan(size(Y));
pred_comb = nan(size(Y));


% Generate predictions from the final integrated regression models.
pred_KL(valid_log)   = predict(mdl1_log, KLlog(valid_log));
pred_eps(valid_log)  = predict(mdl2_log, EPSlog(valid_log));
pred_comb(valid_log) = predict(mdl3_log, Xreg_log(valid_log,:));

% Convert to column vectors.
obs       = Y(:);
pred_KL   = pred_KL(:);
pred_eps  = pred_eps(:);
pred_comb = pred_comb(:);

% Use observations for which all three predictions and the observed
% flux are finite.
% valid = ~(isnan(obs) | ...
%           isnan(pred_KL) | ...
%           isnan(pred_eps) | ...
%           isnan(pred_comb));
% OR
valid_KL = isfinite(Y) & isfinite(pred_KL);
valid_eps = isfinite(Y) & isfinite(pred_eps);
valid_comb = isfinite(Y) & isfinite(pred_comb);

obs       = obs(valid);
pred_KL   = pred_KL(valid);
pred_eps  = pred_eps(valid);
pred_comb = pred_comb(valid);

% Store predictions in a cell array for plotting.
y_pred = {pred_KL, pred_eps, pred_comb};

%% ----------------------------------------------------------------
%  2. CALCULATE R^2 FOR FIGURE 4
% -----------------------------------------------------------------
%
% The displayed R^2 is calculated directly from the observed and
% predicted values.
%
% For an ordinary least-squares regression evaluated on the same
% observations used to fit the model:
%
%       R^2 = corr(Y, Yhat)^2
%
% This is the quantity that should correspond to the Integrated
% row of Table 5.
%
R2_plot = nan(1,3);

for i = 1:3

    r = corr( ...
        obs, ...
        y_pred{i}, ...
        'Rows', 'complete');

    R2_plot(i) = r^2;

end

%% ----------------------------------------------------------------
%  3. OPTIONAL CONSISTENCY CHECK
% -----------------------------------------------------------------
%
% If the final Table 5 integrated-model R^2 values are available in
% the workspace, compare them directly with the Figure 4 values.
%
% The variables R2_KL_int, R2_EPS_int, and R2_Both_int are generated
% by the Table 5 analysis.
%
if exist('R2_KL_int', 'var') && ...
   exist('R2_EPS_int', 'var') && ...
   exist('R2_Both_int', 'var')

    R2_Table5 = [ ...
        R2_KL_int, ...
        R2_EPS_int, ...
        R2_Both_int];

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('TABLE 5 / FIGURE 4 R^2 CONSISTENCY CHECK\n');
    fprintf('============================================================\n');

    fprintf('KL:       Table 5 = %.6f   Figure 4 = %.6f   Difference = %.6e\n', ...
        R2_Table5(1), R2_plot(1), ...
        R2_plot(1) - R2_Table5(1));

    fprintf('epsilon:  Table 5 = %.6f   Figure 4 = %.6f   Difference = %.6e\n', ...
        R2_Table5(2), R2_plot(2), ...
        R2_plot(2) - R2_Table5(2));

    fprintf('KL+eps:   Table 5 = %.6f   Figure 4 = %.6f   Difference = %.6e\n', ...
        R2_Table5(3), R2_plot(3), ...
        R2_plot(3) - R2_Table5(3));

end

%% ----------------------------------------------------------------
%  4. FIGURE 4: PREDICTED vs OBSERVED
% -----------------------------------------------------------------

figure('Color', 'w');

for i = 1:3

    %% -------------------------------------------------------------
    % Density calculation
    % --------------------------------------------------------------

    nbins  = 60;

    xedges = linspace( ...
        min(obs), ...
        max(obs), ...
        nbins);

    yedges = linspace( ...
        min(y_pred{i}), ...
        max(y_pred{i}), ...
        nbins);

    N = histcounts2( ...
        obs, ...
        y_pred{i}, ...
        xedges, ...
        yedges);

    [~, xbin] = histc(obs, xedges);
    [~, ybin] = histc(y_pred{i}, yedges);

    density = nan(size(obs));

    for k = 1:length(obs)

        if xbin(k) > 0 && ...
           ybin(k) > 0 && ...
           xbin(k) < nbins && ...
           ybin(k) < nbins

            density(k) = N(xbin(k), ybin(k));

        end

    end

    %% -------------------------------------------------------------
    % Regression line and prediction interval
    % -------------------------------------------------------------
    %
    % This robust regression is retained for visualization only.
    % Its R^2 is NOT used in the figure title.
    %
    mdl_plot = fitlm( ...
        obs, ...
        y_pred{i}, ...
        'RobustOpts', 'bisquare');

    xfit = linspace( ...
        min(obs), ...
        max(obs), ...
        300)';

    [yfit, yPI] = predict( ...
        mdl_plot, ...
        xfit, ...
        'Prediction', 'observation');

    %% -------------------------------------------------------------
    % Binned mean +/- 1 sigma
    % -------------------------------------------------------------

    nb = 20;

    binEdges = linspace( ...
        min(obs), ...
        max(obs), ...
        nb+1);

    binCenters = 0.5 * ...
        (binEdges(1:end-1) + binEdges(2:end));

    binMean = nan(nb,1);
    binStd  = nan(nb,1);

    for b = 1:nb

        inBin = ...
            obs >= binEdges(b) & ...
            obs <  binEdges(b+1);

        if sum(inBin) > 30

            binMean(b) = ...
                mean(y_pred{i}(inBin));

            binStd(b) = ...
                std(y_pred{i}(inBin));

        end

    end

    %% -------------------------------------------------------------
    % Plot
    % -------------------------------------------------------------

    subplot(1,3,i);

    hold on;
    box on;

    % Density-coloured observations.
    scatter( ...
        obs, ...
        y_pred{i}, ...
        8, ...
        density, ...
        'filled');

    % 1:1 reference line.
    plot( ...
        xfit, ...
        xfit, ...
        'k--', ...
        'LineWidth', 1.4);

    % Robust regression fit retained as a visual trend line.
    plot( ...
        xfit, ...
        yfit, ...
        'r-', ...
        'LineWidth', 2.2);

    % 95% prediction interval.
    fill( ...
        [xfit; flipud(xfit)], ...
        [yPI(:,1); flipud(yPI(:,2))], ...
        [0.5 0.5 0.5], ...
        'FaceAlpha', 0.4, ...
        'EdgeColor', 'none');

    % Binned mean +/- 1 sigma.
    errorbar( ...
        binCenters, ...
        binMean, ...
        binStd, ...
        'ko', ...
        'MarkerFaceColor', 'k', ...
        'LineWidth', 1.2, ...
        'CapSize', 3);

    ylim([0 8]);

    %% -------------------------------------------------------------
    % Axis labels
    % -------------------------------------------------------------

    if i == 2

        xlabel( ...
            'Observed log_{10} Flux [el.cm^{-2}s^{-1}sr^{-1}]', ...
            'FontSize', 13);

    end

    if i == 1

        ylabel( ...
            'Predicted log_{10} Flux [el.cm^{-2}s^{-1}sr^{-1}]', ...
            'FontSize', 13);

    end

    %% -------------------------------------------------------------
    % Title
    % -------------------------------------------------------------
    %
    % IMPORTANT:
    % Use R2_plot, NOT mdl_plot.Rsquared.Ordinary.
    %

    title( ...
        sprintf( ...
            '%s (R^2 = %.2f)', ...
            panelTitle{i}, ...
            R2_plot(i)));

    %% -------------------------------------------------------------
    % Legend
    % -------------------------------------------------------------

    if i == 1

        legend( ...
            {'Data density', ...
             '1:1 Line', ...
             'Regression fit', ...
             'Prediction interval', ...
             'Binned mean \pm1\sigma'}, ...
            'Location', ...
            'northwest');

    end

    %% -------------------------------------------------------------
    % Colorbar
    % -------------------------------------------------------------

    cb = colorbar;

    cb.Label.String = ...
        'Point density (counts per bin)';

    cb.FontSize = 11;

    %% -------------------------------------------------------------
    % Figure formatting
    % -------------------------------------------------------------

    set(gca, 'FontSize', 12);

    grid on;

end

%% ----------------------------------------------------------------
%  3. FIGURE 5: RESIDUAL PLOTS
%% ----------------------------------------------------------------

res_KL   = obs - pred_KL;
res_eps  = obs - pred_eps;
res_comb = obs - pred_comb;

N_obs = length(res_KL);
nbins = 15;

models = {pred_KL,  pred_eps,  pred_comb};
resids = {res_KL,   res_eps,   res_comb};
titles = {'KL', '\epsilon', 'KL+\epsilon'};

figure('Color', 'w');

for i = 1:3

    edges = linspace(min(models{i}), max(models{i}), nbins+1);
    edges = unique(edges);

    [~, ~, bin] = histcounts(models{i}, edges);

    bin_center = nan(nbins, 1);
    res_std    = nan(nbins, 1);
    res_mean   = nan(nbins, 1);
    counts     = nan(nbins, 1);

    for b = 1:nbins
        idx = bin == b;
        if sum(idx) > 5
            bin_center(b) = mean(models{i}(idx));
            res_std(b)    = std(resids{i}(idx));
            res_mean(b)   = mean(resids{i}(idx));
            counts(b)     = sum(idx);
        end
    end

    percent = 100 * counts / N_obs;

    subplot(1, 3, i);
    scatter(models{i}, resids{i}, 10, 'filled', 'MarkerFaceAlpha', 0.2);
    hold on;

    yyaxis left;
    plot(bin_center, res_mean,           'k-',  'LineWidth', 2);
    plot(bin_center, res_mean + res_std, 'r--', 'LineWidth', 1.5);
    plot(bin_center, res_mean - res_std, 'r--', 'LineWidth', 1.5);
    yline(0, 'k:', 'LineWidth', 1.5);
    ylim([-3 5]);

    if i == 1
        ylabel('Residuals');
    end

    yyaxis right;
    bar(bin_center, percent, 0.5, 'FaceAlpha', 0.3);
    ylabel('% of Total Observations');
    ylim([0 max(percent)*1.3]);

    if i == 2
        xlabel('log_{10} Predicted Flux [el.cm^{-2}s^{-1}sr^{-1}]');
    end

    title(titles{i});
    grid on;
    set(gca, 'FontSize', 11);

    if i == 2
        legend({'Residuals', 'Binned Mean', '\pm 1\sigma', '', 'Zero Line', ...
                '% observation'}, 'Location', 'best');
    end

end
