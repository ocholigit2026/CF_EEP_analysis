%% ================================================================
%  load_data.m
%
%  Central data-loading script for the CF–EEP analysis.
%  Run this script once at the start of a MATLAB session before
%  executing any figure or table script.
%
%  ---------------------------------------------------------------
%  OUTPUTS placed in the base workspace
%
%  fluxdata_2009_2019H_4D : [Nvars x Nt x NMLT x NL] particle data
%                            Variable order (row index):
%                              1  lval
%                              2  folon
%                              3  folat
%                              4  mlt
%                              5  pas0
%                              6  pas90
%                              7  mep0e1   (>30  keV, 0°)
%                              8  mep0e2   (>100 keV, 0°)
%                              9  mep0e3   (>300 keV, 0°)
%                             10  mep90e1  (>30  keV, 90°)
%                             11  mep90e2  (>100 keV, 90°)
%                             12  mep90e3  (>300 keV, 90°)
%                             13  PAD
%
%  Flux          : [Nt x NL x NMLT x NE]
%                  Precipitating (0°) electron flux,
%                  re-ordered from the 4D array above.
%
%  cpl_swdata    : timetable – solar wind + coupling function data
%                  Required fields: Ekl, eps, Dst, AE, Bz, v, Kp
%
%  t_hr          : datetime vector (hourly time axis)
%
%  ---------------------------------------------------------------
%  DATA PATHS — update these to match your local directory layout
%
% ================================================================

%% ================================================================
%% FILE PATHS  (edit as needed)
%% ================================================================

PARTICLE_FILE  = fullfile('C:', 'Users', 'mamud', 'DATA HOUSE', ...
                           'poes', 'processed_particle_sw_data_L', ...
                           'fluxdata_2009_2019H.mat');

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
load(PARTICLE_FILE, 'fluxdata_2009_2019H_4D');

fprintf('Loading solar wind data ...\n');
load(SW_FILE, 'swdata_2009_2019H');

fprintf('Loading coupling function data ...\n');
load(COUPLING_FILE, 'coupling_fns_tt');

%% ================================================================
%% COMBINE SOLAR WIND + COUPLING FUNCTIONS
%% ================================================================

cpl_swdata = synchronize(swdata_2009_2019H, coupling_fns_tt);

%% ================================================================
%% EXTRACT AND RESHAPE FLUX ARRAY
%%   Input  : Nvars x Nt x NMLT x NL
%%   Output : Nt x NL x NMLT x NE   (NE = 3 energy channels)
%% ================================================================

Flux = fluxdata_2009_2019H_4D(7:9, :, :, :);   % 0° telescope, E1–E3
Flux = permute(Flux, [2 4 3 1]);                % → Nt x NL x NMLT x NE

%% ================================================================
%% TIME AXIS
%% ================================================================

t_hr = cpl_swdata.Time;

fprintf('Data loaded successfully.\n');
fprintf('  Flux size : %s\n', mat2str(size(Flux)));
fprintf('  Time span : %s  to  %s\n', ...
    datestr(t_hr(1)), datestr(t_hr(end)));
