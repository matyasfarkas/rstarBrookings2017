using DSGE;
using Plots # no need for `using Plots` as that is reexported here

path = dirname(@__FILE__)
horizon  = 80

# Standard model
m10 = Model1010("ss20");
system10 = compute_system(m10)
states_irf10, obs_irf10, pseudo_irf10 = impulse_responses(system10, horizon)
# No FG model
m11 = Model1011("ss20");
system11 = compute_system(m11)
states_irf11, obs_irf11, pseudo_irf11 = impulse_responses(system11, horizon)
# No liquidity and safety model
m12 = Model1012("ss20");
system12 = compute_system(m12)
states_irf12, obs_irf12, pseudo_irf12 = impulse_responses(system12, horizon)

plot(1:horizon,[pseudo_irf10[m10.pseudo_observables[:Forward5YearRealNaturalRate],:, m10.exogenous_shocks[:b_liqp_sh]], pseudo_irf11[ m11.pseudo_observables[:Forward5YearRealNaturalRate],:, m11.exogenous_shocks[:b_liqp_sh]] ],title="IRFs of r* to permanent liquidity shock", label=["Basline model" "No FG model"])
plot!(legend=:bottomright)
savefig( "irf/all/IRF_rstar_to_permanet_liquidity_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[pseudo_irf10[m10.pseudo_observables[:Forward5YearRealNaturalRate],:, m10.exogenous_shocks[:b_safep_sh]], pseudo_irf11[ m11.pseudo_observables[:Forward5YearRealNaturalRate],:, m11.exogenous_shocks[:b_safep_sh]] ],title="IRFs of r*  to permanent safety shock", label=["Basline model" "No FG model"])
plot!(legend=:bottomright)
savefig( "irf/all/IRF_rstar_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf10[ m10.pseudo_observables[:Forward5YearRealNaturalRate],:, m10.exogenous_shocks[:zp_sh ]],pseudo_irf11[ m11.pseudo_observables[:Forward5YearRealNaturalRate],:, m11.exogenous_shocks[:zp_sh ]],pseudo_irf12[ m12.pseudo_observables[:Forward5YearRealNaturalRate],:, m12.exogenous_shocks[:zp_sh ]] ],title="IRFs of r* to permanent technology shock", label=["Basline model" "No FG model" "No convenience yield shocks model"])
plot!(legend=:bottomright)
savefig( "irf/all/IRF_rstar_to_permanet_technology_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf10[ m10.pseudo_observables[:Forward5YearRealNaturalRate],:, m10.exogenous_shocks[:z_sh ]],pseudo_irf11[ m11.pseudo_observables[:Forward5YearRealNaturalRate],:, m11.exogenous_shocks[:z_sh ]],pseudo_irf12[ m12.pseudo_observables[:Forward5YearRealNaturalRate],:, m12.exogenous_shocks[:z_sh ]] ],title="IRFs of r* to transitory technology shock", label=["Basline model" "No FG model" "No convenience yield shocks model"])
plot!(legend=:bottomright)
savefig( "irf/all/IRF_rstar_to_transitory_technology_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf10[ m10.pseudo_observables[:Forward5YearRealNaturalRate],:, m10.exogenous_shocks[:μ_sh ]],pseudo_irf11[ m11.pseudo_observables[:Forward5YearRealNaturalRate],:, m11.exogenous_shocks[:μ_sh ]],pseudo_irf12[ m12.pseudo_observables[:Forward5YearRealNaturalRate],:, m12.exogenous_shocks[:μ_sh ]] ],title="IRFs of r* to investment specific technology shock", label=["Basline model" "No FG model" "No convenience yield shocks model"])
plot!(legend=:bottomright)
savefig( "irf/all/IRF_rstar_to_investment_specific_technology_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf10[ m10.pseudo_observables[:Forward5YearRealNaturalRate],:, m10.exogenous_shocks[:σ_ω_sh ]],pseudo_irf11[ m11.pseudo_observables[:Forward5YearRealNaturalRate],:, m11.exogenous_shocks[:σ_ω_sh ]],pseudo_irf12[ m12.pseudo_observables[:Forward5YearRealNaturalRate],:, m12.exogenous_shocks[:σ_ω_sh ]] ],title="IRFs of r* to risk shock", label=["Basline model" "No FG model" "No convenience yield shocks model"])
plot!(legend=:bottomright)
savefig( "irf/all/IRF_rstar_to_risk_shock.pdf")   # saves the plot from p as a .pdf vector graphic


# policy rate
plot(1:horizon,[obs_irf10[m10.observables[:obs_nominalrate],:, m10.exogenous_shocks[:b_liqp_sh]], obs_irf11[m11.observables[:obs_nominalrate],:, m11.exogenous_shocks[:b_liqp_sh]] ],title="IRFs of the policy rate to permanent liquidity shock", label=["Basline model" "No FG model"])
plot!(legend=:bottomright)
savefig( "irf/all/IRF_robs_to_permanet_liquidity_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[obs_irf10[m10.observables[:obs_nominalrate],:, m10.exogenous_shocks[:b_safep_sh]], obs_irf11[m11.observables[:obs_nominalrate],:, m11.exogenous_shocks[:b_safep_sh]] ],title="IRFs of  the policy rate to permanent safety shock", label=["Basline model" "No FG model"])
plot!(legend=:bottomright)
savefig( "irf/all/IRF_robs_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ obs_irf10[m10.observables[:obs_nominalrate],:, m10.exogenous_shocks[:zp_sh ]],obs_irf11[m11.observables[:obs_nominalrate],:, m11.exogenous_shocks[:zp_sh ]],obs_irf12[m12.observables[:obs_nominalrate],:, m12.exogenous_shocks[:zp_sh ]] ],title="IRFs of  the policy rateto permanent technology shock", label=["Basline model" "No FG model" "No convenience yield shocks model"])
plot!(legend=:bottomright)
savefig( "irf/all/IRF_robs_to_permanet_technology_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ obs_irf10[m10.observables[:obs_nominalrate],:, m10.exogenous_shocks[:z_sh ]],obs_irf11[m11.observables[:obs_nominalrate],:, m11.exogenous_shocks[:z_sh ]],obs_irf12[m12.observables[:obs_nominalrate],:, m12.exogenous_shocks[:z_sh ]] ],title="IRFs of  the policy rateto transitory technology shock", label=["Basline model" "No FG model" "No convenience yield shocks model"])
plot!(legend=:bottomright)
savefig( "irf/all/IRF_robs_to_transitory_technology_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ obs_irf10[m10.observables[:obs_nominalrate],:, m10.exogenous_shocks[:μ_sh ]],obs_irf11[m11.observables[:obs_nominalrate],:, m11.exogenous_shocks[:μ_sh ]],obs_irf12[m12.observables[:obs_nominalrate],:, m12.exogenous_shocks[:μ_sh ]] ],title="IRFs of  the policy rateto investment specific technology shock", label=["Basline model" "No FG model" "No convenience yield shocks model"])
plot!(legend=:bottomright)
savefig( "irf/all/IRF_robs_to_investment_specific_technology_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ obs_irf10[m10.observables[:obs_nominalrate],:, m10.exogenous_shocks[:σ_ω_sh ]],obs_irf11[m11.observables[:obs_nominalrate],:, m11.exogenous_shocks[:σ_ω_sh ]],obs_irf12[m12.observables[:obs_nominalrate],:, m12.exogenous_shocks[:σ_ω_sh ]] ],title="IRFs of  the policy rateto risk shock", label=["Basline model" "No FG model" "No convenience yield shocks model"])
plot!(legend=:bottomright)
savefig( "irf/all/IRF_robs_to_risk_shock.pdf")   # saves the plot from p as a .pdf vector graphic



# obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh]]
# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_safep_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_safep_sh]],obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:b_safep_sh]]  ],title="IRFs to permanent safety shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita"])
# plot!(legend=:bottomright)
# savefig( "irf/m1010/IRF_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic


# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_liqp_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:b_liqp_sh]]   ],title="IRFs to permanent liquidity shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
# plot!(legend=:bottomright)
# savefig( "irf/m1010/IRF_to_permanet_liquidity_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_safep_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_safep_sh]],obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:b_safep_sh]]  ],title="IRFs to permanent safety shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita"])
# plot!(legend=:bottomright)
# savefig( "irf/m1010/IRF_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:zp_sh ]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:zp_sh ]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:zp_sh ]]   ],title="IRFs to permanent technology shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
# plot!(legend=:bottomright)
# savefig( "irf/m1010/IRF_to_permanet_technology_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:π_star_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:π_star_sh]],obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:π_star_sh]]  ],title="IRFs to inflation target shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita"])
# plot!(legend=:bottomright)
# savefig( "irf/m1010/IRF_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:rm_shl6]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:rm_shl6]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:rm_shl6]] ],title="IRFs to 6 period-ahead forward guidance shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate"  "Y:GDP per capita"])
# plot!(legend=:bottomright)
# savefig( "irf/m1010/IRF_to_FG6_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:λ_f_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:λ_f_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:λ_f_sh]]   ],title="IRFs to a price markup shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
# plot!(legend=:bottomright)
# savefig( "irf/m1010/IRF_to_price_markup_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:λ_w_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:λ_w_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:λ_w_sh]]   ],title="IRFs to a wage markup shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
# plot!(legend=:bottomright)
# savefig( "irf/m1010/IRF_to_wage_markup_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# using DSGE;
# path = dirname(@__FILE__)
# m = Model1011("ss20");

# system = compute_system(m)
# horizon =80
# states_irf, obs_irf, pseudo_irf = impulse_responses(system, horizon)

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_liqp_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:b_liqp_sh]]   ],title="IRFs to permanent liquidity shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
# plot!(legend=:bottomright)
# savefig( "irf/m1011/IRF_to_permanet_liquidity_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_safep_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_safep_sh]],obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:b_safep_sh]]  ],title="IRFs to permanent safety shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita"])
# plot!(legend=:bottomright)
# savefig( "irf/m1011/IRF_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:zp_sh ]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:zp_sh ]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:zp_sh ]]   ],title="IRFs to permanent technology shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
# plot!(legend=:bottomright)
# savefig( "irf/m1011/IRF_to_permanet_technology_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:π_star_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:π_star_sh]],obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:π_star_sh]]  ],title="IRFs to inflation target shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita"])
# plot!(legend=:bottomright)
# savefig( "irf/m1011/IRF_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:λ_f_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:λ_f_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:λ_f_sh]]   ],title="IRFs to a price markup shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
# plot!(legend=:bottomright)
# savefig( "irf/m1011/IRF_to_price_markup_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:λ_w_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:λ_w_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:λ_w_sh]]   ],title="IRFs to a wage markup shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
# plot!(legend=:bottomright)
# savefig( "irf/m1011/IRF_to_wage_markup_shock.pdf")   # saves the plot from p as a .pdf vector graphic

using DSGE, Plots;
path = dirname(@__FILE__)
m = Model1012("ss20");

system = compute_system(m)
horizon =80
states_irf, obs_irf, pseudo_irf = impulse_responses(system, horizon)

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_liqp_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_liqp_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:b_liqp_sh]]   ],title="IRFs to permanent liquidity shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
# plot!(legend=:bottomright)
# savefig( "irf/m1012/IRF_to_permanet_liquidity_shock.pdf")   # saves the plot from p as a .pdf vector graphic

# plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:b_safep_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:b_safep_sh]],obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:b_safep_sh]]  ],title="IRFs to permanent safety shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita"])
# plot!(legend=:bottomright)
# savefig( "irf/m1012/IRF_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:zp_sh ]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:zp_sh ]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:zp_sh ]]   ],title="IRFs to permanent technology shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
plot!(legend=:bottomright)
savefig( "irf/m1012/IRF_to_permanet_technology_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:π_star_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:π_star_sh]],obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:π_star_sh]]  ],title="IRFs to inflation target shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita"])
plot!(legend=:bottomright)
savefig( "irf/m1012/IRF_to_permanet_safety_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:λ_f_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:λ_f_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:λ_f_sh]]   ],title="IRFs to a price markup shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
plot!(legend=:bottomright)
savefig( "irf/m1012/IRF_to_price_markup_shock.pdf")   # saves the plot from p as a .pdf vector graphic

plot(1:horizon,[ pseudo_irf[ m.pseudo_observables[:Forward5YearRealNaturalRate],:, m.exogenous_shocks[:λ_w_sh]], obs_irf[ m.observables[:obs_nominalrate],:, m.exogenous_shocks[:λ_w_sh]] , obs_irf[ m.observables[:obs_gdp],:, m.exogenous_shocks[:λ_w_sh]]   ],title="IRFs to a wage markup shock", label=["r*: 5-Year Forward Real Natural Rate" "R: Policy rate" "Y:GDP per capita" ])
plot!(legend=:bottomright)
savefig( "irf/m1012/IRF_to_wage_markup_shock.pdf")   # saves the plot from p as a .pdf vector graphic

