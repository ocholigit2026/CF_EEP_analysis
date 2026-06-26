# CF–EEP Analysis

MATLAB scripts reproducing all figures and tables in the manuscript:

> **[Manuscript title here]**  
> [Author(s)], [Journal], [Year]

---

## Repository structure

```
CF_EEP_analysis/
│
├── load_data.m                          ← Run this first (loads all datasets)
│
├── figures/
│   ├── fig01_MLT_lagged_correlation.m          → Figure 1
│   ├── fig02_global_lagged_correlation.m        → Figure 2
│   ├── fig03_tables1_2_3_integration_window.m  → Figure 3, Tables 1–3
│   ├── fig04_fig05_scatter_residual_plots.m    → Figures 4 & 5
│   └── fig06_AE_mediation_epoch_composite.m    → Figure 6
│
├── tables/
│   └── table4_AE_mediation_analysis.m          → Table 4
│
└── utils/
    ├── bootstrap_CI.m          helper: bootstrap confidence interval
    ├── buildEpochMatrix.m      helper: superposed epoch matrix builder
    └── selectPeakEvents.m      helper: isolate one peak per event cluster
```

---

## Script → output mapping

| Script | Figures | Tables | Notes |
|--------|---------|--------|-------|
| `fig01_MLT_lagged_correlation.m`        | Fig. 1 | — | MLT-sector lagged correlations |
| `fig02_global_lagged_correlation.m`     | Fig. 2 | — | All-MLT global lagged correlations |
| `fig03_tables1_2_3_integration_window.m`| Fig. 3 | Tables 1, 2, 3 | Integration window sweep; instantaneous and peak correlations; regression + VIF |
| `fig04_fig05_scatter_residual_plots.m`  | Figs. 4, 5 | — | Predicted vs. observed; residual diagnostics |
| `table4_AE_mediation_analysis.m`        | — | Table 4 | AE mediation; partial correlations; ΔR² |
| `fig06_AE_mediation_epoch_composite.m`  | Fig. 6 | — | Storm vs. nonstorm superposed epoch composites |

**Table 5** is produced by re-running the regression block in
`fig03_tables1_2_3_integration_window.m` with log-transformed
instantaneous/lagged/integrated coupling inputs and reporting R²
only. No separate script is needed.

---

## How to run

1. Edit the file paths in `load_data.m` to match your local data directory.
2. In MATLAB, run `load_data` to populate the workspace.
3. Add the `utils/` folder to the MATLAB path:
   ```matlab
   addpath('utils')
   ```
4. Run any figure or table script individually, e.g.:
   ```matlab
   run figures/fig01_MLT_lagged_correlation.m
   ```

> All figure scripts depend on variables created by `load_data.m`.  
> `fig04_fig05_scatter_residual_plots.m` additionally depends on the
> regression models (`mdl1`, `mdl2`, `mdl3`) produced by
> `fig03_tables1_2_3_integration_window.m`; run that script first.

---

## Data requirements

The scripts expect three `.mat` files:

| Variable in file | Description |
|-----------------|-------------|
| `fluxdata_2009_2019H_4D` | POES particle flux, 4D array `[Nvars × Nt × NMLT × NL]`, hourly, 2009–2019 |
| `swdata_2009_2019H` | OMNI solar wind timetable, hourly |
| `coupling_fns_tt` | Coupling function timetable (`Ekl`, `eps`, etc.) |

These files are not included in the repository due to size.  
Contact the corresponding author for data access.

---

## Dependencies

- MATLAB R2020b or later (uses `movsum`, `omitmissing`, `tiledlayout`)
- Statistics and Machine Learning Toolbox (`fitlm`, `partialcorr`, `zscore`)
- No additional toolboxes required

---

## License

[Add your licence here, e.g. MIT or CC BY 4.0]
