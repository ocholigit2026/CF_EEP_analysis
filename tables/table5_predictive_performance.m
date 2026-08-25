%% ================================================================
%  TABLE 5: PREDICTIVE PERFORMANCE — INSTANTANEOUS, LAGGED,
%           AND INTEGRATED COUPLING FUNCTIONS
%
%  PURPOSE:
%  Compare R² of regression models using log-transformed coupling
%  functions under three forcing representations:
%    1. Instantaneous
%    2. Predictive-optimal lag (searched over 0–50 h)
%    3. Temporally integrated (peak windows from Table 2)
%
%  ---------------------------------------------------------------
%  REQUIRED WORKSPACE (from load_data.m)
%
%  KL, EPS      : raw coupling functions
%  logJp30      : log10 of >30 keV precipitating flux
%  (logJp100, logJp300 for other channels — change Y below)
%
% ================================================================

%% ================================================================
%% 1. SETTINGS
%% ================================================================

intKL  = 6;     % integration window for KL  (from Table 2)
intEPS = 15;    % integration window for EPS (from Table 2)
maxLag = 50;    % search range for predictive-optimal lag (hours)

offsetKL  = 1e-3;
offsetEPS = 1e-3;

%% ================================================================
%% 2. ALIGN SERIES LENGTH
%% ================================================================

N = min([numel(logJp30), numel(KL), numel(EPS)]);

Y   = logJp30(1:N);
KLv = KL(1:N);
EPSv = EPS(1:N);

%% ================================================================
%% 3. LOG-TRANSFORM RAW COUPLING FUNCTIONS
%% ================================================================

KLlog  = log10(KLv  + offsetKL);
EPSlog = log10(EPSv + offsetEPS);

valid0 = isfinite(Y) & isfinite(KLlog) & isfinite(EPSlog);
fprintf('Valid observations: %d of %d\n', sum(valid0), N);

%% ================================================================
%% 4. INSTANTANEOUS MODELS
%% ================================================================

Y0   = Y(valid0);
KL0  = KLlog(valid0);
EPS0 = EPSlog(valid0);

mdl_KL_inst   = fitlm(KL0,        Y0);
mdl_EPS_inst  = fitlm(EPS0,       Y0);
mdl_Both_inst = fitlm([KL0 EPS0], Y0);

R2_KL_inst   = mdl_KL_inst.Rsquared.Ordinary;
R2_EPS_inst  = mdl_EPS_inst.Rsquared.Ordinary;
R2_Both_inst = mdl_Both_inst.Rsquared.Ordinary;

%% ================================================================
%% 5. PREDICTIVE-OPTIMAL LAG SEARCH
%% ================================================================

R2_KL_lag  = NaN(maxLag+1, 1);
R2_EPS_lag = NaN(maxLag+1, 1);

for lag = 0:maxLag

    if lag == 0
        xKL = KLlog;   xEPS = EPSlog;   YY = Y;
    else
        xKL  = KLlog(1:end-lag);
        xEPS = EPSlog(1:end-lag);
        YY   = Y(1+lag:end);
    end

    vKL  = isfinite(xKL)  & isfinite(YY);
    vEPS = isfinite(xEPS) & isfinite(YY);

    if sum(vKL)  > 2
        m = fitlm(xKL(vKL),   YY(vKL));
        R2_KL_lag(lag+1)  = m.Rsquared.Ordinary;
    end
    if sum(vEPS) > 2
        m = fitlm(xEPS(vEPS), YY(vEPS));
        R2_EPS_lag(lag+1) = m.Rsquared.Ordinary;
    end

end

[~, idxKL]  = max(R2_KL_lag);
[~, idxEPS] = max(R2_EPS_lag);

optLagKL  = idxKL  - 1;
optLagEPS = idxEPS - 1;

fprintf('KL  optimal lag = %d h\n', optLagKL);
fprintf('EPS optimal lag = %d h\n', optLagEPS);

%% ================================================================
%% 6. FIT OPTIMAL-LAG MODELS
%% ================================================================

% KL
X  = KLlog(1:end-optLagKL);
YY = Y(1+optLagKL:end);
v  = isfinite(X) & isfinite(YY);
m  = fitlm(X(v), YY(v));
R2_KL_optlag = m.Rsquared.Ordinary;

% EPS
X  = EPSlog(1:end-optLagEPS);
YY = Y(1+optLagEPS:end);
v  = isfinite(X) & isfinite(YY);
m  = fitlm(X(v), YY(v));
R2_EPS_optlag = m.Rsquared.Ordinary;

% KL + EPS (align to longer lag)
maxOptLag    = max(optLagKL, optLagEPS);
Y_both       = Y(maxOptLag+1:end);
KL_both      = KLlog(maxOptLag-optLagKL+1  : end-optLagKL);
EPS_both     = EPSlog(maxOptLag-optLagEPS+1 : end-optLagEPS);

v = isfinite(Y_both) & isfinite(KL_both) & isfinite(EPS_both);
m = fitlm([KL_both(v) EPS_both(v)], Y_both(v));
R2_Both_optlag = m.Rsquared.Ordinary;

%% ================================================================
%% 7. INTEGRATED MODELS
%% ================================================================

KL_int_log  = log10(movmean(KLv,  [intKL  0], 'omitmissing') + offsetKL);
EPS_int_log = log10(movmean(EPSv, [intEPS 0], 'omitmissing') + offsetEPS);

v = isfinite(Y) & isfinite(KL_int_log);
m = fitlm(KL_int_log(v), Y(v));
R2_KL_int = m.Rsquared.Ordinary;

v = isfinite(Y) & isfinite(EPS_int_log);
m = fitlm(EPS_int_log(v), Y(v));
R2_EPS_int = m.Rsquared.Ordinary;

v = isfinite(Y) & isfinite(KL_int_log) & isfinite(EPS_int_log);
m = fitlm([KL_int_log(v) EPS_int_log(v)], Y(v));
R2_Both_int = m.Rsquared.Ordinary;

%% ================================================================
%% TABLE 5
%% ================================================================

ModelType  = {'Instantaneous'; 'Optimal lag'; 'Integrated'};
R2_KL_col  = [R2_KL_inst;   R2_KL_optlag;   R2_KL_int];
R2_EPS_col = [R2_EPS_inst;  R2_EPS_optlag;  R2_EPS_int];
R2_Both_col = [R2_Both_inst; R2_Both_optlag; R2_Both_int];

Table5 = table(ModelType, R2_KL_col, R2_EPS_col, R2_Both_col, ...
    'VariableNames', {'ModelType', 'KL_R2', 'EPS_R2', 'KL_EPS_R2'});

disp('=== TABLE 5: Predictive Performance (log-transformed inputs) ===');
disp(Table5);

fprintf('\nKL  integration window : %d h\n', intKL);
fprintf('EPS integration window : %d h\n', intEPS);
