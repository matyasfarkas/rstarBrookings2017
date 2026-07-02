# Replication quickstart

This file gives the shortest path to regenerate the paper artifacts. For the
full figure/table map, see `README.md`.

Smoke-tested on June 30, 2026 with Julia 1.5.0. The scripts in the run order
below completed successfully in the local Windows/Parallels setup.

## Environment

Use Julia 1.5.x with the working package environment that includes:

`DSGE`, `Plots`, `StatsPlots`, `CSV`, `DataFrames`, `XLSX`, `HDF5`, `JLD2`,
`FileIO`, `OrderedCollections`, `Measures`, `ClusterManagers`, and
`ModelConstructors`.

Run all commands from the repository root.

## Run order

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

Some charts and tables are composed in Excel from generated CSV/PNG/PDF files.
The main workbooks are:

- `Main results/US/Chart US.xlsx`
- `Main results/US/HLW vs DSGE.xlsx`
- `Main results/EA/Charts EA Final.xlsx`
- `Main results/Raw_Policy_Shocks.xlsx`
- `irf/irf/Paper IRFs.xlsx`

## Smoke test

After running the scripts, these files should exist:

```powershell
$files = @(
  "dsge/Final Paper/Figures/interest_rate_peg_combined_mit_vs_fg_vs_mpplusfg.pdf",
  "dsge/Final Paper/Figures/interest_rate_peg_combined_other_observables.pdf",
  "dsge/Final Paper/Figures/USrstar_wCB_replication.pdf",
  "dsge/Final Paper/Figures/EArstar_wCB_replication.pdf",
  "dsge/Final Paper/Counterfactual/What_if_real_rate_gap_change_HLW_using_policy_rate_shock.pdf",
  "counterfactual/paper/EA_What_if_real_rate_gap_change_HLW_using_policy_rate_shock.pdf",
  "dsge/Final Paper/Figures/Exorbitant privilege 3x2 baseline_plus_delta.pdf",
  "dsge/Final Paper/Figures/EA Case2B no exorbitant privilege 3x2.pdf"
)

$files | ForEach-Object {
  [pscustomobject]@{ File = $_; Exists = Test-Path -LiteralPath $_ }
}
```
