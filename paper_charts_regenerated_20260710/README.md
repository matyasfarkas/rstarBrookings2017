# Regenerated Paper Charts

This folder collects the regenerated chart artifacts for `R-stars Across the
Atlantic - The Role of Policy Expectations`.

Generated on: 2026-07-10

## What Was Regenerated

The documented Julia chart scripts were rerun before collecting these files:

- `dsge/Main_US_rstar_forecasts.jl`
- `dsge/Main_EA_rstar_forecasts.jl`
- `dsge/US_rstar_w_confidence_bands.jl`
- `dsge/EA_rstar_w_confidence_bands.jl`
- `irf/Paper_IRF_Charts_combined_updated.jl`
- `irf/Paper_IRF_Charts_combined_other_observables.jl`
- `counterfactual/Final_US_HLW_RR_gap.jl`
- `counterfactual/Final_EA_HLW_RR_gap.jl`
- `counterfactual/FINAL_Case2b_No_Exorbitant_Privilege.jl`
- `counterfactual/FINAL_Case2d_EA_noCY.jl`

## Folder Layout

- `01_actual_expected_rates_yields_spreads`: Excel workbook sources for Figure 1.
- `02_rstar_and_output_gap`: US and EA r-star/output-gap panel outputs.
- `03_policy_innovations`: Excel/data sources for the policy-innovation chart.
- `04_policy_shock_contributions`: HVD workbook/image sources for policy-shock contribution charts.
- `05_interest_rate_peg_policy_shocks`: monetary-policy peg IRF figure.
- `06_interest_rate_peg_convenience_yield_shocks`: convenience-yield shock peg IRFs.
- `07_other_observables_and_spreads`: policy/spread/OIS/convenience-yield IRFs.
- `08_uncertainty_bands_and_starfish`: confidence-band and starfish charts.
- `09_HLW_counterfactuals`: US and EA HLW counterfactual outputs.
- `10_convenience_yield_estimates`: convenience-yield estimates and related chart sources.
- `11_US_exorbitant_privilege_counterfactual`: US exorbitant-privilege counterfactual outputs.
- `appendix_observables`: observable-variable workbook/image sources.
- `appendix_rstar_decompositions`: r-star decomposition images and workbook sources.
- `rate_domain_variant_not_standard_paper`: optional CY interest-rate-domain outputs from the separate rate-domain exercise.

## Alignment Note

PDF and PNG chart files were copied directly from the regenerated script outputs.
They were not resized, cropped, or recompressed, so panel alignment and layout are
preserved exactly as produced by the paper scripts.

Excel-composed figures are included as workbook sources because those charts are
stored inside workbooks rather than regenerated as standalone PDFs by the Julia
pipeline.

## Verification

The collection step copied 83 files and reported 0 missing files. See:

- `_sources/copied_files_manifest.csv`
- `_sources/missing_files.csv`
