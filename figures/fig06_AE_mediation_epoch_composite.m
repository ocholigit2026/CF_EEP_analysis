%% ================================================================
%  FIGURE 6: SUPERPOSED EPOCH ANALYSIS - STORM vs NONSTORM
%
%  6-panel stacked figure:
%    Panel 1 - KL reconnection electric field [mV/m]
%    Panel 2 - IMF Bz [nT]
%    Panel 3 - Solar wind speed v [km/s]
%    Panel 4 - AE index [nT]
%    Panel 5 - Dst index [nT]
%    Panel 6 - Log10 precipitating flux (>30 keV)
%
%  Storm    = min(Dst) <= -30 nT within +/- 12 h of KL peak
%  Nonstorm = otherwise
%  Bootstrap 95% confidence bands (n = 1000).
%
%  ---------------------------------------------------------------
%  REQUIRED WORKSPACE (from load_data.m)
%
%  cpl_swdata : timetable with fields Ekl, eps, Dst, AE, Bz, v, flux30_0deg
%
%  Helper functions (in utils/):
%    bootstrap_CI.m
%    buildEpochMatrix.m
%    selectPeakEvents.m
%    classifyStormNonstorm.m
%
% ================================================================

%% ----------------------------------------------------------------
%  1. SETTINGS
%% ----------------------------------------------------------------

epochWin  = 96;
threshPct = 95;
nBoot     = 1000;
alpha     = 0.05;
minSep    = 60;
stormDst  = -30;
stormWin  = 12;

%% ----------------------------------------------------------------
%  2. EXTRACT VARIABLES
%% ----------------------------------------------------------------

KL    = cpl_swdata.Ekl;
EPS   = cpl_swdata.eps;
Dst   = cpl_swdata.Dst;
FluxP = cpl_swdata.flux30_0deg;

%% ----------------------------------------------------------------
%  3. EVENT SELECTION
%% ----------------------------------------------------------------

KLthr  = prctile(KL,  threshPct);
EPSthr = prctile(EPS, threshPct);

KL_events  = selectPeakEvents(find(KL  >= KLthr),  KL,  minSep);
EPS_events = selectPeakEvents(find(EPS >= EPSthr), EPS, minSep);

fprintf('KL events:  %d\n', length(KL_events));
fprintf('EPS events: %d\n', length(EPS_events));

%% ----------------------------------------------------------------
%  4. STORM / NONSTORM CLASSIFICATION
%% ----------------------------------------------------------------

[KL_storm_idx,  KL_nonstorm_idx]  = classifyStormNonstorm(KL_events,  Dst, stormDst, stormWin);
[EPS_storm_idx, EPS_nonstorm_idx] = classifyStormNonstorm(EPS_events, Dst, stormDst, stormWin);

fprintf('Storm KL events:     %d\n', length(KL_storm_idx));
fprintf('Nonstorm KL events:  %d\n', length(KL_nonstorm_idx));
fprintf('Storm EPS events:    %d\n', length(EPS_storm_idx));
fprintf('Nonstorm EPS events: %d\n', length(EPS_nonstorm_idx));

%% ----------------------------------------------------------------
%  5. BUILD EPOCH MATRICES
%% ----------------------------------------------------------------

cmp     = @(evts, ser) buildEpochMatrix(evts, ser, epochWin);
trimNaN = @(M) M(~any(isnan(M), 2), :);

P_KL_storm    = trimNaN(cmp(KL_storm_idx,    FluxP)');
P_KL_nonstorm = trimNaN(cmp(KL_nonstorm_idx, FluxP)');

Bz_storm_KL  = trimNaN(cmp(KL_storm_idx,    cpl_swdata.Bz)');
Bz_quiet_KL  = trimNaN(cmp(KL_nonstorm_idx, cpl_swdata.Bz)');

v_storm_KL   = trimNaN(cmp(KL_storm_idx,    cpl_swdata.v)');
v_quiet_KL   = trimNaN(cmp(KL_nonstorm_idx, cpl_swdata.v)');

AE_storm_KL  = trimNaN(cmp(KL_storm_idx,    cpl_swdata.AE)');
AE_quiet_KL  = trimNaN(cmp(KL_nonstorm_idx, cpl_swdata.AE)');

Dst_storm_KL = trimNaN(cmp(KL_storm_idx,    cpl_swdata.Dst)');
Dst_quiet_KL = trimNaN(cmp(KL_nonstorm_idx, cpl_swdata.Dst)');

KL_storm_mat = trimNaN(cmp(KL_storm_idx,    cpl_swdata.Ekl)');
KL_quiet_mat = trimNaN(cmp(KL_nonstorm_idx, cpl_swdata.Ekl)');

%% ----------------------------------------------------------------
%  6. BOOTSTRAP CONFIDENCE INTERVALS
%% ----------------------------------------------------------------

bootCI = @(X) bootstrap_CI(X, nBoot, alpha);

[mean_storm_KL,  st_lo_KL,    st_hi_KL]       = bootCI(KL_storm_mat);
[mean_quiet_KL,  qt_lo_KL,    qt_hi_KL]       = bootCI(KL_quiet_mat);

[mean_bz_storm_KL, bz_st_lo_KL, bz_st_hi_KL] = bootCI(Bz_storm_KL);
[mean_bz_quiet_KL, bz_qt_lo_KL, bz_qt_hi_KL] = bootCI(Bz_quiet_KL);

[mean_v_storm_KL,  v_st_lo_KL,  v_st_hi_KL]  = bootCI(v_storm_KL);
[mean_v_quiet_KL,  v_qt_lo_KL,  v_qt_hi_KL]  = bootCI(v_quiet_KL);

[mean_ae_storm_KL, ae_st_lo_KL, ae_st_hi_KL] = bootCI(AE_storm_KL);
[mean_ae_quiet_KL, ae_qt_lo_KL, ae_qt_hi_KL] = bootCI(AE_quiet_KL);

[mean_dst_storm_KL, dst_st_lo_KL, dst_st_hi_KL] = bootCI(Dst_storm_KL);
[mean_dst_quiet_KL, dst_qt_lo_KL, dst_qt_hi_KL] = bootCI(Dst_quiet_KL);

[mean_stKL, ~, ~] = bootCI(P_KL_storm);
[mean_qtKL, ~, ~] = bootCI(P_KL_nonstorm);

%% ----------------------------------------------------------------
%  7. PLOTTING
%% ----------------------------------------------------------------

eW = (-epochWin:epochWin);

shade = @(ax, eW, lo, hi, col, alph) fill(ax, ...
    [eW fliplr(eW)], [lo fliplr(hi)], col, ...
    'EdgeColor', 'none', 'FaceAlpha', alph);

figure('Color', 'w');
tiledlayout(6, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Panel 1: KL
nexttile; hold on;
shade(gca, eW, st_lo_KL,  st_hi_KL,  'r', 0.10);
shade(gca, eW, qt_lo_KL,  qt_hi_KL,  'b', 0.50);
plot(eW, mean_storm_KL, 'r', 'LineWidth', 2);
plot(eW, mean_quiet_KL, 'b', 'LineWidth', 2);
xline(0, 'k--');
legend('', '', 'Storm', 'Nonstorm', 'Location', 'northeast', 'Orientation', 'horizontal');
ylabel('E_{kl} [mV/m]');
xticklabels(''); grid on; box on; xlim([-48 eW(end)]);

% Panel 2: Bz
nexttile; hold on;
shade(gca, eW, bz_st_lo_KL, bz_st_hi_KL, 'r', 0.20);
shade(gca, eW, bz_qt_lo_KL, bz_qt_hi_KL, 'b', 0.20);
plot(eW, mean_bz_storm_KL, 'r', 'LineWidth', 3);
plot(eW, mean_bz_quiet_KL, 'b', 'LineWidth', 3);
xline(0, 'k--');
ylabel('B_z [nT]');
xticklabels(''); grid on; box on; xlim([-48 eW(end)]);

% Panel 3: Solar wind speed
nexttile; hold on;
shade(gca, eW, v_st_lo_KL, v_st_hi_KL, 'r', 0.20);
shade(gca, eW, v_qt_lo_KL, v_qt_hi_KL, 'b', 0.20);
plot(eW, mean_v_storm_KL, 'r', 'LineWidth', 3);
plot(eW, mean_v_quiet_KL, 'b', 'LineWidth', 3);
xline(0, 'k--');
ylabel('v [km/s]');
xticklabels(''); grid on; box on; xlim([-48 eW(end)]); ylim([350 550]);

% Panel 4: AE
nexttile; hold on;
shade(gca, eW, ae_st_lo_KL, ae_st_hi_KL, 'r', 0.20);
shade(gca, eW, ae_qt_lo_KL, ae_qt_hi_KL, 'b', 0.20);
plot(eW, mean_ae_storm_KL, 'r', 'LineWidth', 3);
plot(eW, mean_ae_quiet_KL, 'b', 'LineWidth', 3);
xline(0, 'k--');
ylabel('AE [nT]');
xticklabels(''); grid on; box on; xlim([-48 eW(end)]);

% Panel 5: Dst
nexttile; hold on;
shade(gca, eW, dst_st_lo_KL, dst_st_hi_KL, 'r', 0.20);
shade(gca, eW, dst_qt_lo_KL, dst_qt_hi_KL, 'b', 0.20);
plot(eW, mean_dst_storm_KL, 'r', 'LineWidth', 3);
plot(eW, mean_dst_quiet_KL, 'b', 'LineWidth', 3);
xline(0, 'k--');
ylabel('Dst [nT]');
xticklabels(''); grid on; box on; xlim([-48 eW(end)]);

% Panel 6: Precipitating flux
nexttile; hold on;
plot(eW, log10(mean_stKL), 'r', 'LineWidth', 2);
plot(eW, log10(mean_qtKL), 'b', 'LineWidth', 2);
xline(0, 'k--');
ylabel('Log_{10}J_p [el.cm^{-2}s^{-1}sr^{-1}]');
xlabel('Epoch Time [Hours]');
grid on; box on; xlim([-48 eW(end)]);
