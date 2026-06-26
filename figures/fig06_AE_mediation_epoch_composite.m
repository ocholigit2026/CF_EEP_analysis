%% ================================================================
%  FIGURE 6: SUPERPOSED EPOCH ANALYSIS — STORM vs NONSTORM
%            COMPOSITE RESPONSE TO COUPLING FUNCTION PEAKS
%
%  PURPOSE:
%  Provide a qualitative view of the AE mediation role by comparing
%  the composite time-evolution of key variables during KL-peak
%  events classified as geomagnetic storm vs. nonstorm periods.
%
%  OUTPUT:
%  6-panel stacked figure (tiledlayout):
%    Panel 1 – Kan-Lee coupling function  E_kl
%    Panel 2 – IMF Bz
%    Panel 3 – Solar wind speed v
%    Panel 4 – AE index
%    Panel 5 – Dst index
%    Panel 6 – Log10 precipitating flux (>30 keV, 0° telescope)
%
%  Storm = min(Dst) ≤ −30 nT within ±12 h of KL peak
%  Nonstorm = otherwise
%
%  Bootstrap 95% confidence bands are shown.
%
%  ---------------------------------------------------------------
%  REQUIRED WORKSPACE VARIABLES (loaded via load_data.m)
%
%  cpl_swdata : timetable with fields
%               Ekl, eps, Dst, AE, Bz, v, Kp,
%               flux30_0deg, flux30_90deg
%
% ================================================================

%% ================================================================
%% USER SETTINGS
%% ================================================================

epochWin  = 96;    % ± hours for epoch window
threshPct = 95;    % top-percentile threshold for event selection
nBoot     = 2000;  % bootstrap iterations
alpha     = 0.05;  % 95% CI

minSep    = 60;    % minimum separation between events (hours)
stormDst  = -30;   % Dst threshold for storm classification (nT)
stormWin  = 12;    % window around onset for Dst check (± hours)

%% ================================================================
%% EXTRACT VARIABLES
%% ================================================================

KL    = cpl_swdata.Ekl;
EPS   = cpl_swdata.eps;
Dst   = cpl_swdata.Dst;
FluxP = cpl_swdata.flux30_0deg;
FluxT = cpl_swdata.flux30_90deg;
t     = cpl_swdata.Time;

epochWindow = (-epochWin:epochWin)';   % epoch axis (hours)
ntau        = length(epochWindow);

%% ================================================================
%% EVENT SELECTION (Top-percentile, isolated KL peaks)
%% ================================================================

KLthr  = prctile(KL,  threshPct);
EPSthr = prctile(EPS, threshPct);

KL_candidates  = find(KL  >= KLthr);
EPS_candidates = find(EPS >= EPSthr);

KL_events  = selectPeakEvents(KL_candidates,  KL,  minSep);
EPS_events = selectPeakEvents(EPS_candidates, EPS, minSep);

fprintf('KL events:  %d\n', length(KL_events));
fprintf('EPS events: %d\n', length(EPS_events));

%% ================================================================
%% STORM / NONSTORM CLASSIFICATION
%% ================================================================

KL_storm_idx    = [];
KL_nonstorm_idx = [];

for i = 1:length(KL_events)
    idx    = KL_events(i);
    wStart = max(1, idx - stormWin);
    wEnd   = min(length(Dst), idx + stormWin);

    if min(Dst(wStart:wEnd)) <= stormDst
        KL_storm_idx(end+1) = idx;
    else
        KL_nonstorm_idx(end+1) = idx;
    end
end

EPS_storm_idx    = [];
EPS_nonstorm_idx = [];

for i = 1:length(EPS_events)
    idx    = EPS_events(i);
    wStart = max(1, idx - stormWin);
    wEnd   = min(length(Dst), idx + stormWin);

    if min(Dst(wStart:wEnd)) <= stormDst
        EPS_storm_idx(end+1) = idx;
    else
        EPS_nonstorm_idx(end+1) = idx;
    end
end

fprintf('Storm KL events:    %d\n', length(KL_storm_idx));
fprintf('Nonstorm KL events: %d\n', length(KL_nonstorm_idx));
fprintf('Storm EPS events:   %d\n', length(EPS_storm_idx));
fprintf('Nonstorm EPS events:%d\n', length(EPS_nonstorm_idx));

%% ================================================================
%% BUILD EPOCH MATRICES
%% ================================================================

cmp = @(evts, ser) buildEpochMatrix(evts, ser, epochWin);
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

KL_storm_KL  = trimNaN(cmp(KL_storm_idx,    cpl_swdata.Ekl)');
KL_quiet_KL  = trimNaN(cmp(KL_nonstorm_idx, cpl_swdata.Ekl)');

%% ================================================================
%% BOOTSTRAP CONFIDENCE INTERVALS
%% ================================================================

bootCI = @(X) bootstrap_CI(X, nBoot, alpha);

[mean_storm_KL,  st_lo_KL,    st_hi_KL]    = bootCI(KL_storm_KL);
[mean_quiet_KL,  qt_lo_KL,    qt_hi_KL]    = bootCI(KL_quiet_KL);

[mean_bz_storm_KL,  bz_st_lo_KL, bz_st_hi_KL] = bootCI(Bz_storm_KL);
[mean_bz_quiet_KL,  bz_qt_lo_KL, bz_qt_hi_KL] = bootCI(Bz_quiet_KL);

[mean_v_storm_KL,   v_st_lo_KL,  v_st_hi_KL]  = bootCI(v_storm_KL);
[mean_v_quiet_KL,   v_qt_lo_KL,  v_qt_hi_KL]  = bootCI(v_quiet_KL);

[mean_ae_storm_KL,  ae_st_lo_KL, ae_st_hi_KL] = bootCI(AE_storm_KL);
[mean_ae_quiet_KL,  ae_qt_lo_KL, ae_qt_hi_KL] = bootCI(AE_quiet_KL);

[mean_dst_storm_KL, dst_st_lo_KL, dst_st_hi_KL] = bootCI(Dst_storm_KL);
[mean_dst_quiet_KL, dst_qt_lo_KL, dst_qt_hi_KL] = bootCI(Dst_quiet_KL);

[mean_stKL, lo_stKL, hi_stKL] = bootCI(P_KL_storm);
[mean_qtKL, lo_qtKL, hi_qtKL] = bootCI(P_KL_nonstorm);

%% ================================================================
%% PLOTTING (6-panel stacked layout)
%% ================================================================

eW = epochWindow';   % row vector for fill() calls

figure('Color', 'w')
tiledlayout(6, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Helper: shaded confidence band
shadeBand = @(ax, eW, lo, hi, col, alpha_v) ...
    fill(ax, [eW fliplr(eW)], [lo fliplr(hi)], col, ...
         'EdgeColor', 'none', 'FaceAlpha', alpha_v);

%% --- Panel 1: KL ---
nexttile; hold on;
shadeBand(gca, eW, st_lo_KL,  st_hi_KL,  'r', 0.10);
shadeBand(gca, eW, qt_lo_KL,  qt_hi_KL,  'b', 0.50);
plot(eW, mean_storm_KL, 'r', 'LineWidth', 2);
plot(eW, mean_quiet_KL, 'b', 'LineWidth', 2);
xline(0, 'k--');
legend('', '', 'Storm', 'Nonstorm', 'Location', 'northeast', 'Orientation', 'horizontal');
ylabel('E_{kl} [mV/m]');
xticklabels(''); grid on; box on;
xlim([-48 eW(end)])

%% --- Panel 2: Bz ---
nexttile; hold on;
shadeBand(gca, eW, bz_st_lo_KL, bz_st_hi_KL, 'r', 0.20);
shadeBand(gca, eW, bz_qt_lo_KL, bz_qt_hi_KL, 'b', 0.20);
plot(eW, mean_bz_storm_KL, 'r', 'LineWidth', 3);
plot(eW, mean_bz_quiet_KL, 'b', 'LineWidth', 3);
xline(0, 'k--');
ylabel('B_z [nT]');
xticklabels(''); grid on; box on;
xlim([-48 eW(end)])

%% --- Panel 3: Solar wind speed ---
nexttile; hold on;
shadeBand(gca, eW, v_st_lo_KL, v_st_hi_KL, 'r', 0.20);
shadeBand(gca, eW, v_qt_lo_KL, v_qt_hi_KL, 'b', 0.20);
plot(eW, mean_v_storm_KL, 'r', 'LineWidth', 3);
plot(eW, mean_v_quiet_KL, 'b', 'LineWidth', 3);
xline(0, 'k--');
ylabel('v [km/s]');
xticklabels(''); grid on; box on;
xlim([-48 eW(end)]); ylim([350 550])

%% --- Panel 4: AE ---
nexttile; hold on;
shadeBand(gca, eW, ae_st_lo_KL, ae_st_hi_KL, 'r', 0.20);
shadeBand(gca, eW, ae_qt_lo_KL, ae_qt_hi_KL, 'b', 0.20);
plot(eW, mean_ae_storm_KL, 'r', 'LineWidth', 3);
plot(eW, mean_ae_quiet_KL, 'b', 'LineWidth', 3);
xline(0, 'k--');
ylabel('AE [nT]');
xticklabels(''); grid on; box on;
xlim([-48 eW(end)])

%% --- Panel 5: Dst ---
nexttile; hold on;
shadeBand(gca, eW, dst_st_lo_KL, dst_st_hi_KL, 'r', 0.20);
shadeBand(gca, eW, dst_qt_lo_KL, dst_qt_hi_KL, 'b', 0.20);
plot(eW, mean_dst_storm_KL, 'r', 'LineWidth', 3);
plot(eW, mean_dst_quiet_KL, 'b', 'LineWidth', 3);
xline(0, 'k--');
ylabel('Dst [nT]');
xticklabels(''); grid on; box on;
xlim([-48 eW(end)])

%% --- Panel 6: Precipitating flux ---
nexttile; hold on;
plot(eW, log10(mean_stKL), 'r', 'LineWidth', 2);
plot(eW, log10(mean_qtKL), 'b', 'LineWidth', 2);
xline(0, 'k--');
ylabel('Log_{10}J_p [el.cm^{-2}s^{-1}sr^{-1}]');
xlabel('Epoch Time [Hours]');
grid on; box on;
xlim([-48 eW(end)])
