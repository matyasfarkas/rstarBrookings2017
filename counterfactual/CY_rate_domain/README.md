# Convenience Yield Rate-Domain Counterfactuals

This folder contains isolated copies for revisiting the US exorbitant privilege
counterfactual with the convenience-yield wedge expressed in the short-rate
domain.

The conversion used throughout is:

```julia
CY_rate_domain_APR = 4 * ( -1 / Xi_r ) * CY
```

where `Xi_r` is the policy-rate coefficient in the Euler equation:

```julia
Xi_r = (1 - h * exp(-z_star)) / (sigma_c * (1 + h * exp(-z_star)))
```

## Scripts

- `Extract_ConvenienceYield_PosteriorMode_rate_domain.jl` extracts the US and EA
  posterior-mode convenience-yield components and exports them in rate-domain
  APR units.
- `FINAL_Case2b_No_Exorbitant_Privilege_rate_domain.jl` is the copied Case2b
  counterfactual script. It leaves the original script untouched and redirects
  all outputs to `dsge/Final Paper/Figures/CY_rate_domain`.

## Suggested Run Order

```powershell
julia counterfactual/CY_rate_domain/Extract_ConvenienceYield_PosteriorMode_rate_domain.jl
julia counterfactual/CY_rate_domain/FINAL_Case2b_No_Exorbitant_Privilege_rate_domain.jl
```
