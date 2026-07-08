# R-stars Across the Atlantic replication package

This repository is the replication package for *R-stars Across the Atlantic -
The Role of Policy Expectations* by Matyas Farkas, Zoltan Jakab, and Jesper
Linde, IMF Working Paper, June 2026.

The package contains the Julia scripts, model inputs, chart workbooks, and final
paper outputs used to reproduce the United States and euro-area results in the
paper. It is intentionally kept close to the research working directory so that
the published figures, intermediate outputs, and historical scripts remain
traceable.

## Core relationship to the original r-star paper

This project builds directly on the FRBNY/Brookings `rstarBrookings2017`
replication code for *Safety, Liquidity, and the Natural Rate of Interest* by
Marco Del Negro, Domenico Giannone, Marc Giannoni, and Andrea Tambalotti,
*Brookings Papers on Economic Activity*, Spring 2017.

The original paper and code base provide the DSGE infrastructure, r-star
measurement framework, safe/liquid asset convenience-yield block, and many of
the model conventions used here. The Atlantic paper extends that framework to
study policy expectations, US-euro-area comparisons, interest-rate peg
counterfactuals, and convenience-yield exercises.

Original r-star paper:
[*Safety, Liquidity, and the Natural Rate of Interest*](https://www.brookings.edu/bpea-articles/safety-liquidity-and-the-natural-rate-of-interest/)

Original FRBNY repository:
[`FRBNY-DSGE/rstarBrookings2017`](https://github.com/FRBNY-DSGE/rstarBrookings2017)

## What this package reproduces

- US and euro-area DSGE estimates of r-star and related observables.
- Confidence-band charts and post-pandemic r-star projections.
- Interest-rate peg impulse responses under surprise, anticipated, and mixed
  monetary-policy shocks.
- Alternative policy-peg charts for OIS, AAA spreads, BAA spreads, and
  convenience-yield responses.
- Counterfactual paths under HLW real-rate gaps.
- US and euro-area convenience-yield counterfactuals, including the
  no-exorbitant-privilege exercises.
- Excel-composed charts and tables used in the final paper package.

## Recommended branch

Use the `replication` branch for the paper-specific replication guide, script
map, smoke-test checklist, and small runtime fixes needed in the local
Julia 1.5.0 environment.

```powershell
git checkout replication
```

The `master` branch keeps the fork visible as a descendant of the original
FRBNY repository, while pointing readers to the Atlantic-paper replication
workflow.

## Software

The final scripts were smoke-tested with Julia 1.5.0 on Windows/Parallels paths
pointing to a Mac-mounted working tree.

Julia packages used by the final scripts include:

- `DSGE`
- `Plots`, `StatsPlots`, `Measures`
- `CSV`, `DataFrames`, `XLSX`
- `HDF5`, `JLD2`, `FileIO`
- `OrderedCollections`, `ClusterManagers`, `ModelConstructors`
- Julia standard libraries: `Dates`, `LinearAlgebra`, `Statistics`

Excel is required for paper charts and tables that were composed from model
outputs in workbooks. There is no root `Project.toml`; use a Julia environment
with the packages above installed.

## Quick replication run order

From the repository root on the `replication` branch, run:

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

For a shorter command-only guide, see `REPLICATION.md` on the `replication`
branch.

## Main directories

- `dsge/`: DSGE model scripts, inputs, forecast outputs, and final paper
  figures.
- `irf/`: impulse-response and policy-peg chart scripts.
- `counterfactual/`: HLW, convenience-yield, and no-exorbitant-privilege
  counterfactuals.
- `Main results/`: Excel workbooks, raw chart inputs, and final chart
  workbooks for the US and euro area.
- `tvar/` and `plot/`: legacy material inherited from the original
  `rstarBrookings2017` replication package.

## Key final outputs

Representative outputs are written to:

- `dsge/Final Paper/Figures`
- `dsge/Final Paper/Counterfactual`
- `counterfactual/paper`
- `Main results/US`
- `Main results/EA`
- `irf/irf/Paper IRFs.xlsx`

Running the scripts may overwrite existing generated outputs in those
directories. If exact archival preservation matters, copy the repository before
rerunning the full sequence.

## Citation

If you use this replication package, please cite:

Farkas, Matyas, Zoltan Jakab, and Jesper Linde. 2026. *R-stars Across the
Atlantic - The Role of Policy Expectations*. IMF Working Paper.

Please also cite the foundational r-star paper and code base:

Del Negro, Marco, Domenico Giannone, Marc Giannoni, and Andrea Tambalotti.
2017. "Safety, Liquidity, and the Natural Rate of Interest." *Brookings Papers
on Economic Activity*, Spring 2017, 235-294.

## License and disclaimer

This repository is a research fork of the FRBNY `rstarBrookings2017` code base.
The original license terms in `LICENSE` continue to apply. The original code
was provided by the Federal Reserve Bank of New York on an "as is" basis,
without warranties or conditions of any kind.
