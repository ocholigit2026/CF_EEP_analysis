%% ================================================================
%  load_data.m
%
%  Central data-loading script for the CF–EEP analysis.
%  Run this script ONCE at the start of a MATLAB session before
%  executing any figure or table script.
%
%  ---------------------------------------------------------------
%  WORKSPACE OUTPUTS
%
%  fluxdata_2009_2019H_4D_new  [nVar x nTime x nMLT x nL]
%
%    Variable index (row):
%      7   mep0e1  (>30  keV, 0°)
%      8   mep0e2  (>100 keV, 0°)
%      9   mep0e3  (>300 keV, 0°)
%     10   mep90e1 (>30  keV, 90°)
%     11   mep90e2 (>100 keV, 90°)
%     12   mep90e3 (>300 keV, 90°)
%
%    Dimension ordering:
%      nVar = 12 (selected variables)
%      nTime = number of hourly steps (2009–2019)
%      nMLT  = 24  (1-h MLT bins, 0–23)
%      nL    = 57  (0.25 L bins, 1.00–15.00)
%
%  ParticleHourly   timetable: Jp30, Jp100, Jp300,
%                              logJp30, logJp100, logJp300
%  ParticleMLT30/100/300  timetables: time x 24 MLT columns
%
%  cpl_swdata       timetable: Ekl, eps, Dst, AE, Bz, v, Kp, ...
%  t_hr             datetime vector (hourly)
%  Lidx             L = 4–6 index vector (13:21)
%
%  ---------------------------------------------------------------
%  NOTE ON DATA FILE VERSIONS
%
%  This script loads fluxdata_2009_2019H_4D_new produced by
%  meped_preprocessing.m.  The OLD file fluxdata_2009_2019H.mat
%  uses a different dimension ordering (variables x time x L x MLT)
%  and is NO LONGER USED by any analysis script in this repository.
%
%  ---------------------------------------------------------------
%  FILE PATHS  (edit to match your local directory layout)
%
% ================================================================

%% ================================================================
%% ADD UTILS TO PATH
%% ================================================================

addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

%% ================================================================
%% FILE PATHS
%% ================================================================

PARTICLE_FILE  = fullfile('C:', 'Users', 'mamud', 'DATA HOUSE', ...
                           'poes', 'raw_data', ...
                           'MEPED_particle_reproduced_optimized.mat');

SW_FILE        = fullfile('C:', 'Users', 'mamud', 'DATA HOUSE', ...
                           'OMNIWEB DATA', 'hourly', ...
                           'swdataTT_2009_2019H_filtered.mat');

COUPLING_FILE  = fullfile('C:', 'Users', 'mamud', 'DATA HOUSE', ...
                           'energy_couple_fns', ...
                           'energy_couple_data_2009_2019.mat');

%% ================================================================
%% LOAD
%% ================================================================

fprintf('Loading particle data ...\n');
load(PARTICLE_FILE, ...
    'fluxdata_2009_2019H_4D_new', ...
    'ParticleHourly', ...
    'ParticleMLT30', 'ParticleMLT100', 'ParticleMLT300', ...
    'L_edges', 'MLT_edges', 'Lidx', 't_hr');

fprintf('Loading solar wind data ...\n');
load(SW_FILE, 'swdata_2009_2019H');

fprintf('Loading coupling function data ...\n');
load(COUPLING_FILE, 'coupling_fns_tt');

%% ================================================================
%% COMBINE SOLAR WIND + COUPLING FUNCTIONS
%% ================================================================

cpl_swdata = synchronize(swdata_2009_2019H, coupling_fns_tt);

%% ================================================================
%% CONVENIENCE VARIABLES
%% ================================================================

% L-shell and MLT grids (consistent with preprocessing)
L_values   = 1:0.25:15;     % 57 values
MLT_values = 0:23;          % 24 values

% L = 4–6 index range (Lidx = 13:21)
% Already loaded from PARTICLE_FILE; recomputed here for transparency
Lidx = find(L_values >= 4 & L_values <= 6);

% Global 1-D flux series (log-transformed, from ParticleHourly)
logJp30  = ParticleHourly.logJp30;
logJp100 = ParticleHourly.logJp100;
logJp300 = ParticleHourly.logJp300;

% Add global flux series to cpl_swdata for convenience
cpl_swdata.flux30_0deg  = ParticleHourly.Jp30;
cpl_swdata.flux100_0deg = ParticleHourly.Jp100;
cpl_swdata.flux300_0deg = ParticleHourly.Jp300;

% Raw coupling proxies
KL  = cpl_swdata.Ekl(:);
EPS = cpl_swdata.eps(:);

% Standardised coupling proxies (used by fig01, fig02, fig03)
KL_z  = (KL  - mean(KL,  'omitnan')) ./ std(KL,  'omitnan');
EPS_z = (EPS - mean(EPS, 'omitnan')) ./ std(EPS, 'omitnan');

%% ================================================================
%% VALIDATION REPORT
%% ================================================================

[nVar, nTime, nMLT, nL] = size(fluxdata_2009_2019H_4D_new);

fprintf('\nData loaded successfully.\n');
fprintf('  4-D array : [%d x %d x %d x %d]  (var x time x MLT x L)\n', ...
    nVar, nTime, nMLT, nL);
fprintf('  Time span : %s  to  %s\n', ...
    datestr(t_hr(1)), datestr(t_hr(end)));
fprintf('  Jp30  finite: %d / %d\n', ...
    sum(isfinite(ParticleHourly.Jp30)), nTime);
fprintf('  Jp100 finite: %d / %d\n', ...
    sum(isfinite(ParticleHourly.Jp100)), nTime);
fprintf('  Jp300 finite: %d / %d\n', ...
    sum(isfinite(ParticleHourly.Jp300)), nTime);
