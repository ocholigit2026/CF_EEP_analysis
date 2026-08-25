# CF–EEP Analysis

MATLAB analysis scripts for the manuscript:

> **Solar Wind–Magnetosphere Coupling Functions as Proxies for Energetic Electron Precipitation: Evidence for Time-Scale Dependence and Internal Modulation**
> Ocholi, M., Akala, A.O., Oyeyemi, E.O., & Reeves, G.D.
> *Journal of Geophysical Research: Space Physics*, 2026
> Code archived at: https://doi.org/10.5281/zenodo.20915979

---

## Repository structure

```
CF_EEP_analysis/
│
├── load_data.m                          ← Run once per session before analyses
│
├── figures/
│   ├── fig01_global_lagged_correlation.m       → Figure 1
│   ├── fig02_MLT_lagged_correlation.m          → Figure 2
│   ├── fig03_tables1_2_3_integration_window.m → Figure 3, Tables 1–3
│   ├── fig04_fig05_scatter_residual_plots.m   → Figures 4 & 5
│   └── fig06_AE_mediation_epoch_composite.m   → Figure 6
│
├── tables/
│   ├── table4_AE_mediation_analysis.m  → Table 4
│   └── table5_predictive_performance.m → Table 5
│
├── utils/
│   ├── bootstrap_CI.m          helper: bootstrap confidence interval
│   ├── buildEpochMatrix.m      helper: superposed epoch matrix builder
│   ├── selectPeakEvents.m      helper: one peak per event cluster
│   └── surrogate_null_test.m   → circular-shift surrogate null analysis
│
└── preprocessing/
    └── meped_preprocessing.m   ← Standalone: reproduces particle data file
```

---

## How to run

```matlab
% 1. Edit file paths in load_data.m to match your data directory

% 2. Load all data into workspace (run once per session)
%    This also adds utils/ to the MATLAB path automatically.
run load_data.m

% 4. Run any figure or table script, e.g.:
run figures/fig01_global_lagged_correlation.m
run tables/table4_AE_mediation_analysis.m
run utils/surrogate_null_test.m
```

**Dependency order within the analysis flow:**
- `load_data.m` must run before any figure or table script
- `fig03_tables1_2_3_integration_window.m` must run before `fig04_fig05_scatter_residual_plots.m`
  (fig04/05 use `mdl1`, `mdl2`, `mdl3`, `KL_inth`, `EPS_inth`, `Xreg` from fig03)
- All other scripts are independent of each other

---

## Preprocessing (standalone)

`preprocessing/meped_preprocessing.m` documents how raw multi-satellite MEPED
observations were aggregated into the hourly L-shell × MLT gridded data structure.
It is provided for transparency and reproducibility and is **not** part of the
main analysis flow.

Run it only if you need to reproduce the processed particle data file
(`MEPED_particle_reproduced_optimized.mat`) from the raw satellite timetables.
`load_data.m` loads this pre-existing file directly.

---

## Script → output mapping

| Script | Figures | Tables |
|--------|---------|--------|
| `fig01_global_lagged_correlation.m` | Fig. 1 | — |
| `fig02_MLT_lagged_correlation.m` | Fig. 2 | — |
| `fig03_tables1_2_3_integration_window.m` | Fig. 3 | Tables 1, 2, 3 |
| `fig04_fig05_scatter_residual_plots.m` | Figs. 4, 5 | — |
| `table4_AE_mediation_analysis.m` | — | Table 4 |
| `table5_predictive_performance.m` | — | Table 5 |
| `fig06_AE_mediation_epoch_composite.m` | Fig. 6 | — |
| `surrogate_null_test.m` | Null figures | — |

---

## Workspace variables provided by load_data.m

| Variable | Description |
|----------|-------------|
| `fluxdata_2009_2019H_4D_new` | `[nVar × nTime × nMLT × nL]` particle flux array |
| `ParticleHourly` | Timetable: `Jp30`, `Jp100`, `Jp300`, `logJp30`, `logJp100`, `logJp300` |
| `ParticleMLT30/100/300` | Timetables: time × 24 MLT columns |
| `cpl_swdata` | Synchronised solar wind + coupling function timetable |
| `t_hr` | Hourly datetime vector |
| `KL`, `EPS` | Raw coupling functions |
| `KL_z`, `EPS_z` | Standardised coupling functions |
| `logJp30`, `logJp100`, `logJp300` | Log10 precipitating flux series |
| `L_values`, `MLT_values` | Spatial grids |
| `Lidx` | L = 4–6 index vector (13:21) |

---

## Data requirements

| File | Contents | Source |
|------|----------|--------|
| `MEPED_particle_reproduced_optimized.mat` | Processed particle data | Produced by `meped_preprocessing.m` |
| `swdataTT_2009_2019H_filtered.mat` | Hourly solar wind timetable | NASA OMNI/SPDF |
| `energy_couple_data_2009_2019.mat` | Coupling function timetable | Derived |

Raw data files are not included due to size. Processed datasets are
available from the corresponding author upon reasonable request.

---

## Dependencies

- MATLAB R2020b or later
- Statistics and Machine Learning Toolbox (`fitlm`, `partialcorr`, `zscore`, `corr`)

---

## License

Creative Commons Attribution 4.0 International (CC BY 4.0)
