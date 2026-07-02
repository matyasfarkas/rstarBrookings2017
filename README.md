# R-stars Across the Atlantic replication package

This branch is a replication-package guide for the IMF Working Paper
`R-stars Across the Atlantic - The Role of Policy Expectations` by Matyas
Farkas, Zoltan Jakab, and Jesper Linde, June 2026.

The repository began as the FRBNY/Brookings `rstarBrookings2017` code base and
has been extended with United States and euro-area model runs, impulse-response
figures, counterfactuals, Excel chart workbooks, and final paper outputs. This
README documents the scripts and workbooks needed to reproduce the paper
artifacts without reorganizing or deleting the existing research files.

For the shortest command-only guide, see `REPLICATION.md`. The quickstart was
smoke-tested on June 30, 2026 with Julia 1.5.0 in the local Windows/Parallels
environment; all listed scripts completed successfully after the small runtime
path/scope fixes included on this branch.

## What is in scope

The replication package preserves the existing repository contents. The
intended reproducible path is:

- Run the final Julia scripts listed below to regenerate model-based CSV, PDF,
  PNG, and XLSX outputs.
- Use the Excel workbooks listed below for charts/tables that were composed in
  Excel from the model outputs.
- Treat `dsge/Final Paper/Figures`, `dsge/Final Paper/Counterfactual`,
  `counterfactual/paper`, and `Main results` as the final paper artifact
  directories.

Running the scripts may overwrite existing output files in those directories.
If exact archival preservation matters, copy the repository before rerunning the
full sequence.

## Software

The scripts were checked on Julia 1.5.0 on Windows/Parallels paths pointing to a
Mac-mounted working tree.

Julia packages used by the final scripts include:

- `DSGE`
- `Plots`, `StatsPlots`, `Measures`
- `CSV`, `DataFrames`, `XLSX`
- `HDF5`, `JLD2`, `FileIO`
- `OrderedCollections`, `ClusterManagers`, `ModelConstructors`
- Julia standard libraries: `Dates`, `LinearAlgebra`, `Statistics`

Excel is required for the workbook-composed charts and tables. There is no root
`Project.toml` in this package; use the working Julia environment that contains
the packages above.

## Key inputs

The final scripts rely on precomputed input data, parameter modes, and forecast
objects already stored in the repository.

Important inputs include:

- `dsge/input_data`
- `dsge/output_data/m1010/ss20/estimate/raw/paramsmode_vint=250825.h5`
- `dsge/output_data/m1010/ss24/estimate/raw/paramsmode_vint=250115.h5`
- `dsge/output_data/m1010/ss24/estimate/raw/paramsmode_vint=250116.h5`
- `Main results/US/Ex_post_real_rate_gaps.csv`
- `Main results/EA/Ex_post_real_rate_gaps.csv`
- `counterfactual/rstarobs_mode`
- `counterfactual/rstarobs_plus_fg`
- `counterfactual/rstarobs_plus_convyieldobs`

## Recommended run order

From the repository root, run the final scripts with Julia:

```powershell
julia "dsge/Main_US_rstar_forecasts.jl"
julia "dsge/Main_EA_rstar_forecasts.jl"
julia "dsge/US_rstar_w_confidence_bands.jl"
julia "dsge/EA_rstar_w_confidence_bands.jl"

julia "irf/Paper_IRF_Charts_combined_updated.jl"
julia "irf/Paper_IRF_Charts_combined_other_observables.jl"

julia "counterfactual/Final_US_HLW_RR_gap.jl"
julia "counterfactual/Final_EA_HLW_RR_gap.jl"
julia "counterfactual/FINAL_Case2b_No_Exorbitant_Privilege.jl"
julia "counterfactual/FINAL_Case2d_EA_noCY.jl"
```

The main estimation scripts under `dsge/Final Paper/Main_full*.jl` and
`dsge/Final Paper/EA_Alternative_YC*.jl` are the historical-decomposition and
forecast-generation scripts. They are heavier than the plotting/counterfactual
scripts and may take materially longer.

## Paper artifact map

| Paper item | Reproduction source | Main outputs |
| --- | --- | --- |
| Figure 1. Actual and expected policy rates, long-term yields, and spreads | Excel-composed from `Main results/US/Chart US.xlsx`, `Main results/EA/Charts EA Final.xlsx`, and raw/input data under `dsge/input_data` and `Main results/*/raw data` | Excel chart exports in `Main results/US` and `Main results/EA` |
| Figure 2. Model estimates of r* and output gaps | `dsge/Main_US_rstar_forecasts.jl`, `dsge/Main_EA_rstar_forecasts.jl`; final panel assembly is stored in `dsge/Final Paper/Figures` | `US_rstar_final_chart.png`, `EA_rstar_final_chart.png`, `US_rstar_final_chart_2016Q1_panel.*`, `EA_rstar_final_chart_2016Q1_panel.*` |
| Table 1. Drivers of filtered r* | Historical-decomposition outputs from `dsge/Final Paper/Main_full*.jl`, `dsge/Final Paper/EA_Alternative_YC*.jl`, and Excel aggregation | `Main results/US/HVD/*.csv`, `Main results/EA/todelete/*shock_decomposition.csv`, `Main results/HVD_diff.xlsx`, `Main results/EA/todelete/EA_HVD_diff.xlsx` |
| Figure 3. Estimated policy innovations and contribution to policy rates | `dsge/Final Paper/Main_full.jl`; Excel chart source `Main results/Raw_Policy_Shocks.xlsx` | `dsge/Final Paper/m1010_histstd_monetary_shocks.csv`, chart workbook output |
| Figure 4. Contributions of monetary-policy shocks to inflation and output | Historical-decomposition CSVs and Excel chart workbooks | `Main results/US/HVD/*.csv`, `Main results/EA/todelete/*shock_decomposition.csv`, chart workbook output |
| Figure 5. Interest-rate peg with contemporaneous and anticipated monetary shocks | `irf/Paper_IRF_Charts_combined_updated.jl` | `dsge/Final Paper/Figures/interest_rate_peg_combined_mit_vs_fg_vs_mpplusfg.pdf`, `irf/irf/Paper IRFs.xlsx` |
| Figure 6. Interest-rate peg with convenience-yield shocks | `irf/Paper_IRF_Charts_combined_updated.jl` | `IRF_rate_peg_with_permanent_liquidity_shock_period0.pdf`, `IRF_rate_peg_with_permanent_and_transitory_liquidity_shocks_period0.pdf`, `IRF_rate_peg_with_permanent_and_transitory_safety_shocks_period0.pdf` |
| Figure 7. Impact of interest-rate peg on interest rates and spreads | `irf/Paper_IRF_Charts_combined_other_observables.jl` | `interest_rate_peg_combined_other_observables.*`, `interest_rate_peg_liquidity_other_observables.*`, `interest_rate_peg_combined_convenience_yield.*`, `interest_rate_peg_liquidity_convenience_yield.*` |
| Figure 8. Uncertainty bands and projected post-pandemic r* estimates | `dsge/US_rstar_w_confidence_bands.jl`, `dsge/EA_rstar_w_confidence_bands.jl`, plus `dsge/Main_US_rstar_forecasts.jl` and `dsge/Main_EA_rstar_forecasts.jl` | `USrstar_wCB_replication.*`, `EArstar_wCB_replication.*`, `US_rstar_starfish_chart*.png`, `EA_rstar_starfish_chart*.png` |
| Figure 9. Counterfactual simulation under HLW r* since 2020Q4 | `counterfactual/Final_US_HLW_RR_gap.jl`, `counterfactual/Final_EA_HLW_RR_gap.jl` | US outputs in `dsge/Final Paper/Counterfactual`; EA outputs in `counterfactual/paper` |
| Figure 10. US and EA convenience-yield estimates | `counterfactual/FINAL_Case2b_No_Exorbitant_Privilege.jl`, `counterfactual/FINAL_Case2d_EA_noCY.jl`; Excel/chart assembly may use the exported CSV | `ConvenienceYield_fullsample_US_EA.csv`, `EA Case2B no exorbitant privilege 3x2.*` |
| Figure 11. Counterfactual with EA convenience-yield innovations in the United States | `counterfactual/FINAL_Case2b_No_Exorbitant_Privilege.jl` | `Exorbitant privilege 3x2 baseline_plus_delta.*`, `US no exorbitant privilege 3x2 baseline_plus_delta.*`, related CSV exports |
| Appendix Figures II.1-II.2. Observable variables | Excel/input-data charts | `Main results/US/Chart US.xlsx`, `Main results/EA/Charts EA Final.xlsx` |
| Appendix Figures IV.1-IV.2. r* decompositions | Excel-composed from HLW and DSGE decomposition outputs | `Main results/US/HLW vs DSGE.xlsx`, `Main results/EA/Charts EA Final.xlsx`, `Main results/US/*decomposition*.png`, `Main results/EA/*decomposition*.png` |

## Excel workbooks

Some final paper graphics were composed in Excel rather than directly saved by
Julia. The main workbook sources are:

- `Main results/US/Chart US.xlsx`
- `Main results/US/HLW vs DSGE.xlsx`
- `Main results/US/Ex_post_real_rate_gaps.xlsx`
- `Main results/EA/Charts EA Final.xlsx`
- `Main results/Rstar_w_withoutFG.xlsx`
- `Main results/Data_correlation_with_policy_rate.xlsx`
- `Main results/Raw_Policy_Shocks.xlsx`
- `irf/irf/Paper IRFs.xlsx`

Temporary Office lock files such as `~$*.xlsx` and LibreOffice lock files are
ignored by `.gitignore`; existing tracked lock files are legacy artifacts.

## Validation checklist

After running the scripts, verify that these representative outputs exist:

```powershell
$files = @(
  "dsge/Final Paper/Figures/interest_rate_peg_combined_mit_vs_fg_vs_mpplusfg.pdf",
  "dsge/Final Paper/Figures/IRF_rate_peg_with_permanent_liquidity_shock_period0.pdf",
  "dsge/Final Paper/Figures/interest_rate_peg_combined_other_observables.pdf",
  "dsge/Final Paper/Figures/interest_rate_peg_liquidity_other_observables.pdf",
  "dsge/Final Paper/Figures/USrstar_wCB_replication.pdf",
  "dsge/Final Paper/Figures/EArstar_wCB_replication.pdf",
  "dsge/Final Paper/Counterfactual/What_if_real_rate_gap_change_HLW_using_policy_rate_shock.pdf",
  "counterfactual/paper/EA_What_if_real_rate_gap_change_HLW_using_policy_rate_shock.pdf",
  "dsge/Final Paper/Figures/Exorbitant privilege 3x2 baseline_plus_delta.pdf",
  "dsge/Final Paper/Figures/EA Case2B no exorbitant privilege 3x2.pdf",
  "Main results/US/Chart US.xlsx",
  "Main results/EA/Charts EA Final.xlsx"
)

$files | ForEach-Object {
  [pscustomobject]@{ File = $_; Exists = Test-Path -LiteralPath $_ }
}
```

All entries should report `Exists = True`.

## Notes for maintainers

- This branch is documentation-oriented. It does not remove legacy files,
  intermediate outputs, or exploratory scripts.
- `.gitignore` is conservative: it ignores new transient files and generated
  cache directories, but it does not untrack outputs already committed to the
  repository.
- Several paths contain spaces, especially under `dsge/Final Paper`; quote paths
  when running scripts from a shell.
- Some scripts contain absolute Windows/Mac-mounted paths from the original
  research environment. If a script fails on another machine, first update those
  paths locally rather than changing the analytical logic.

## License and original code base

The original repository included replication files for `Safety, Liquidity, and
the Natural Rate of Interest` by Marco Del Negro, Domenico Giannone, Marc
Giannoni, and Andrea Tambalotti, Brookings Papers on Economic Activity, Spring
2017. The license terms in `LICENSE` continue to apply.
